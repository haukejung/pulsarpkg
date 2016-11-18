def check_object_type(obj, obj_type, allowNone=False):
    """
    Raises an exception if an object is not of a certain type
    :param obj: The object to test
    :param obj_type: Type the object should have
    """
    if not allowNone and obj is None:
        raise TypeError('{0} is None, when it shouldn\'t be'.format(obj))
    else:
        if obj is not None and not isinstance(obj, obj_type):
            raise TypeError('{0} is not a "{1}" object'.format(obj, obj_type))