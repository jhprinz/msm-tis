import logging
from uuid import UUID
from weakref import WeakValueDictionary

from openpathsampling.mongodb.base import StorableObject
from openpathsampling.mongodb.cache import MaxCache, Cache, NoCache, \
    WeakLRUCache
from openpathsampling.mongodb.proxy import LoaderProxy

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class HashedList(dict):
    def __init__(self):
        super(HashedList, self).__init__()
        dict.__init__(self)
        self._list = []

    def append(self, key):
        dict.__setitem__(self, key, len(self))
        self._list.append(key)

    # noinspection PyCallByClass
    def extend(self, t):
        l = len(self)
        map(lambda x, y: dict.__setitem__(self, x, y), t, range(l, l + len(t)))
        self._list.extend(t)

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        self._list[value] = key

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def index(self, key):
        return self._list[key]

    def mark(self, key):
        if key not in self:
            dict.__setitem__(self, key, -2)

    def unmark(self, key):
        if key in self:
            dict.__delitem__(self, key)

    def clear(self):
        dict.clear(self)
        self._list = []

    @property
    def list(self):
        return self._list


class ObjectStore(StorableObject):
    """
    Base Class for storing complex objects in a netCDF4 file. It holds a
    reference to the store file.`

    Attributes
    ----------
    content_class : :obj:`openpathsampling.netcdfplus.base.StorableObject`
        a reference to the class type to be stored using this Storage. Must be
        subclassed from :obj:`openpathsampling.netcdfplus.base.StorableObject`
    cache : :py:class:`openpathsampling.netcdfplus.cache.Cache`
        a dictionary that holds references to all stored elements by index
        or string for named objects. This is only used for cached access
        if caching is not `False`. Must be of type
        :obj:`openpathsampling.netcdfplus.base.StorableObject` or subclassed.

    """
    _restore_non_initial_attr = False

    allowed_types = [
        'int', 'float', 'long', 'str', 'bool',
        'numpy.float32', 'numpy.float64',
        'numpy.int8', 'numpy.inf16', 'numpy.int32', 'numpy.int64',
        'numpy.uint8', 'numpy.uinf16', 'numpy.uint32', 'numpy.uint64',
        'index', 'length', 'uuid'
    ]

    default_store_chunk_size = 256

    class DictDelegator(object):
        def __init__(self, store, dct):
            self.name = store.name + '_'
            self.dct = dct

        def __getitem__(self, item):
            return self.dct[self.name + item]

        def __contains__(self, item):
            return (self.name + item) in self.dct

    def name_delegate(self, dct):
        return ObjectStore.DictDelegator(self, dct)

    default_cache = 10000

    def __init__(self, name, content_class):
        """

        Parameters
        ----------
        name : str
        content_class : class

        Notes
        -----
        Usually you want caching, but limited. Recommended is to use an LRUCache
        with a reasonable maximum number of objects that depends on the typical
        number of objects to cache and their size

        The class that takes care of storing data in a file is called a
        `Storage`, so the netCDF+ subclassed `Storage` is a storage.
        The classes that know how to load and save an object from the storage
        are called `Store`, like ObjectStore, SampleStore, etc...

        The difference between `json` and `jsonobj` is subtle. Consider
        storing a complex object. Then there are two ways to do that.
        1. `json`: Store a reference to the object (provided) it is stored and
        2. `jsonobj`: serialize the object and only use references for contained
        objects. All inner objects will always be stored using references.
        The only exception is using nestable. Consider objects that contain
        references to objects of the same type, like e.g. operations in an
        equation (2*3 + 3). Each operation represents a value but each
        operation needs values to operate on. To save such an object you have
        again two options:
        1. `nestable=False`. Store all single objects and always reference
        the contained objects. For an equation that would mean to store several
        objects `op1 = plus(op2, 3), op2 = times(2, 3)`. Since this is correct
        though not intuitive you can also use
        2. `nestable=True`. Store all the serialized objects nested into one
        object (string). For our example this corresponds to
        `plus(times(2,3), 3)`.

        """

        super(ObjectStore, self).__init__()
        self._storage = None
        self.content_class = content_class
        self.cache = NoCache()
        self._free = set()
        self._cached_all = False
        self._created = False

        self.name = name

        self.attribute_list = {}
        self.cv = {}

        # This will not be stored since its information is contained in the
        # dimension names
        self._dimension_name_store = None

        self.variables = dict()
        self.units = dict()

        self.index = None

        self.proxy_index = WeakValueDictionary()

        if self.content_class is not None \
                and not issubclass(self.content_class, StorableObject):
            raise ValueError(
                'Content class "%s" must be subclassed from StorableObject.' %
                self.content_class.__name__)

        self.fallback_store = None

    def is_created(self):
        return self._created

    def to_dict(self):
        return {
            'content_class': self.content_class,
            'name': self.name
        }

    def register_fallback(self, store):
        self.fallback_store = store

    def register(self, storage):
        """
        Associate the object store to a specific storage with a given name

        Parameters
        ----------
        storage : :class:`openpathsampling.netcdfplus.NetCDFPlus`
            the storage to be associated with

        """
        self._storage = storage
        self.name = self.name

        self.index = self.create_uuid_index()
        self._document = storage.db[self.name]

    def create_uuid_index(self):
        return HashedList()

    def restore(self):
        self.load_indices()

    def load_indices(self):
        self.index.clear()
        self.index.extend(
            [int(UUID(x)) for x in self._document.distinct('_id')])

    @property
    def storage(self):
        """Return the associated storage object

        Returns
        -------

        :class:`openpathsampling.netcdfplus.NetCDFPlus`
            the referenced storage object
        """

        if self._storage is None:
            raise RuntimeError(
                'A storage needs to be added to this store to be used! '
                'Use .register() to do so.')

        return self._storage

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'store.%s[%s] : %s' % (
            self.name,
            self.content_class.__name__ if self.content_class is not None else
            'None/ANY',
            str(len(self)) + ' object(s)'
        )

    @property
    def simplifier(self):
        """
        Return the simplifier instance used to create JSON serialization

        Returns
        -------
        :class:`openpathsampling.netcdfplus.dictify.StorableObjectJSON`
            the simplifier object used in the associated storage

        """
        return self.storage.simplifier

    def set_caching(self, caching):
        """
        Set the caching mode for this store

        Parameters
        ----------
        caching : :class:`openpathsampling.netcdfplus.Cache`

        """
        if caching is None:
            caching = self.default_cache

        if caching is True:
            caching = MaxCache()
        elif caching is False:
            caching = NoCache()
        elif type(caching) is int:
            caching = WeakLRUCache(caching)

        if isinstance(caching, Cache):
            self.cache = caching.transfer(self.cache)

    def idx(self, obj):
        """
        Return the index in this store for a given object

        Parameters
        ----------
        obj : :class:`openpathsampling.netcdfplus.base.StorableObject`
            the object that can be stored in this store for which its index is
            to be returned

        Returns
        -------
        int or `None`
            The integer index of the given object or `None` if it is not
            stored yet
        """
        return self.index[obj.__uuid__]

    def __iter__(self):
        """
        Add iteration over all elements in the storage
        """
        # we want to iterator in the order object were saved!
        for uuid in self.index._list:
            yield self.load(uuid)

    def __len__(self):
        """
        Return the number of stored objects

        Returns
        -------
        int
            number of stored objects

        """
        return self._document.count()

    # def write(self, variable, idx, obj, attribute=None):
    #     if attribute is None:
    #         attribute = variable
    #
    #     var = self.vars[variable]
    #     val = getattr(obj, attribute)
    #
    #     var[int(idx)] = val
    #
    #     if var.var_type.startswith('lazy'):
    #         proxy = var.store.proxy(val)
    #         if isinstance(obj, LoaderProxy):
    #             # for a loader proxy apply it to the real object
    #             setattr(obj.__subject__, attribute, proxy)
    #         else:
    #             setattr(obj, attribute, proxy)

    def proxy(self, item):
        """
        Return a proxy of a object for this store

        Parameters
        ----------
        item : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            or int The item or index that points to an object in this store
            and to which a proxy is requested.

        Returns
        -------

        """
        if item is None:
            return None

        tt = type(item)
        if tt is long:
            idx = item
        elif tt in [str, unicode]:
            if item[0] == '-':
                return None
            idx = int(UUID(item))
        else:
            idx = item.__uuid__

        return LoaderProxy(self, idx)

    def __contains__(self, item):
        if item.__uuid__ in self.index:
            return True

        if self.fallback_store is not None and item in self.fallback_store:
            return True

        if self.storage.fallback is not None and item in self.storage.fallback:
            return True

        return False

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        try:
            if type(item) is int:
                if item < 0:
                    item += len(self)
                return self.load(item)
            elif type(item) is str or type(item) is long:
                return self.load(item)
            elif type(item) is slice:
                return [self.load(idx)
                        for idx in range(*item.indices(len(self)))]
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return iter(self)
        except KeyError:
            return None

    def get(self, item):
        try:
            return self[item]
        except KeyError:
            return None

    def _load(self, idx):
        obj = self.storage.simplifier.from_simple_dict(
            self._document.find_one({'_id': str(UUID(int=idx))}))
        return obj

    def clear_cache(self):
        """Clear the cache and force reloading"""

        self.cache.clear()
        self._cached_all = False

    def cache_all(self):
        """Load all samples as fast as possible into the cache"""
        if not self._cached_all:
            idxs = range(len(self))
            jsons = self.variables['json'][:]

            [self.add_single_to_cache(i, j) for i, j in zip(
                idxs,
                jsons)]

            self._cached_all = True

    def _save(self, obj, idx):
        dct = self.storage.simplifier.to_simple_dict(obj)
        self._document.insert(dct)

    @property
    def last(self):
        """
        Returns the last generated trajectory. Useful to continue a run.

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the last stored object in this store
        """
        return self.load(len(self) - 1)

    @property
    def first(self):
        """
        Returns the first stored object.

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the actual first stored object
        """
        return self.load(0)

    def free(self):
        """
        Return the number of the next free index for this store

        Returns
        -------
        index : int
            the number of the next free index in the storage.
            Used to store a new object.
        """

        idx = len(self)

        return idx

    def initialize(self):
        """
        Initialize the associated storage to allow for object storage. Mainly
        creates an index dimension with the name of the object.
        """

        self._created = True

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the object to be loaded

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the loaded object
        """

        if type(idx) is long:
            if idx in self.index:
                n_idx = self.index[idx]
            else:
                if self.fallback_store is not None:
                    return self.fallback_store.load(idx)
                elif self.storage.fallback is not None:
                    return self.storage.fallback.stores[self.name].load(idx)
                else:
                    raise ValueError(
                        'str %s not found in storage or fallback' % idx)

        elif type(idx) is not int:
            raise ValueError((
                'indices of type "%s" are not allowed in named storage '
                '(only str and int)') % type(idx).__name__
            )
        else:
            n_idx = int(idx)

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            logger.debug('Found IDX #' + str(idx) + ' in cache. Not loading!')
            return obj

        except KeyError:
            pass

        logger.debug(
            'Calling load object of type `%s` @ IDX #%d' %
            (self.content_class.__name__, n_idx))

        if n_idx >= len(self):
            logger.warning(
                'Trying to load from IDX #%d > number of object %d' %
                (n_idx, len(self)))
            return None
        elif n_idx < 0:
            logger.warning((
                'Trying to load negative IDX #%d < 0. '
                'This should never happen!!!') % n_idx)
            raise RuntimeError(
                'Loading of negative int should result in no object. '
                'This should never happen!')
        else:
            obj = self._load(idx)

        logger.debug(
            'Calling load object of type %s and IDX # %d ... DONE' %
            (self.content_class.__name__, n_idx))

        if obj is not None:
            # update cache there might have been a change due to naming
            self.cache[n_idx] = obj

            logger.debug(
                'Try loading UUID object of type %s and IDX # %d ... DONE' %
                (self.content_class.__name__, n_idx))

        logger.debug(
            'Finished load object of type %s and IDX # %d ... DONE' %
            (self.content_class.__name__, n_idx))

        return obj

    @staticmethod
    def reference(obj):
        return obj.__uuid__

    def remember(self, obj):
        """
        Tell a store that an obj should be assumed as stored

        This is useful, if you do not want to store an object in a specific
        store. Especially to make sure attributes are not stored multiple times

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be fake stored

        """
        self.index.mark(obj.__uuid__)

    def forget(self, obj):
        """
        This will revert remembering non-stored objects.

        Stored objects cannot be forgotten

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be forgotten

        """

        self.index.unmark(obj.__uuid__)

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """
        uuid = obj.__uuid__

        if uuid in self.index:
            # has been saved so quit and do nothing
            if not self.index[uuid] == -1:
                return self.reference(obj)

            # numbers other than -1 are reserved for other things

        if isinstance(obj, LoaderProxy):
            if obj._store is self:
                # is a proxy of a saved object so do nothing
                return uuid
            else:
                # it is stored but not in this store so we try storing the
                # full attribute which might be still in cache or memory
                # if that is not the case it will be stored again. This can
                # happen when you load from one store save to another. And load
                # again after some time while the cache has been changed and try
                # to save again the loaded object. We will not explicitly store
                # a table that matches objects between different storages.
                return self.save(obj.__subject__)

        if self.fallback_store is not None and \
                self.storage.exclude_from_fallback:
            if obj in self.fallback_store:
                return self.reference(obj)

        elif self.storage.fallback is not None and \
                self.storage.exclude_from_fallback:
            if obj in self.storage.fallback:
                return self.reference(obj)

        if not isinstance(obj, self.content_class):
            raise ValueError((
                'This store can only store object of base type "%s". Given '
                'obj is of type "%s". You might need to use another store.')
                % (self.content_class, obj.__class__.__name__)
            )

        n_idx = len(self.index)

        # mark as saved so circular dependencies will not cause infinite loops
        self.index.append(uuid)

        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(n_idx))

        try:
            self._save(obj, n_idx)
            self.cache[n_idx] = obj

        except:
            # in case we did not succeed remove the mark as being saved
            del self.index[uuid]
            raise

        return self.reference(obj)

    def __setitem__(self, key, value):
        """
        Enable saving using __setitem__

        """
        self.save(value, key)

    def add_single_to_cache(self, idx, json):
        """
        Add a single object to cache by json

        Parameters
        ----------
        idx : int
            the index where the object was stored
        json : str
            json string the represents a serialized version of the stored object
        """

        if idx not in self.cache:
            obj = self.simplifier.from_json(json)

            # self._get_id(idx, obj)

            self.cache[idx] = obj
            self.index[obj.__uuid__] = idx

            return obj
