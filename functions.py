def check_object_type(obj, obj_type):
    """
    Raises an exception if an object is not of a certain type
    :param obj: The object to test
    :param obj_type: Type the object should have
    """
    if not isinstance(obj, obj_type):
            TypeError('{0} is not a "{1}" object'.format(obj, obj_type))