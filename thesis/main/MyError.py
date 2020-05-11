class PostProcessError(Exception):
    pass


class SubProcessTimeout(PostProcessError):

    def __init__(self, timeout: float) -> None:
        self.message = "Workers took longer than {t} seconds to finish".format(t=timeout)
        return None


class DataframeEmptyError(PostProcessError):
    def __init__(self, df_name: str) -> None:
        self.message = "Dataframe {name} was empty".format(name=df_name)
        return None


class ParameterError(Exception):
    pass


class DuplicateParameterError(ParameterError):

    def __init__(self, parameter_name: str) -> None:
        self.message = "Multiple options for physical parameter {name} found in parameter set." \
                       "Check consistency of parameter definition.".format(name=parameter_name)
        return None


class DuplicateCollectionError(ParameterError):
    def __init__(self, collection_name: str) -> None:
        self.message = "Multiple options for collection {name} found in parameter set." \
                       "Check consistency of parameter definition.".format(name=collection_name)


class CollectionNotFoundInParameterSet(ParameterError):

    def __init__(self, collection_name: str) -> None:
        self.message = "The Collection {name} was not found." \
                       "Check consistency of parameter definition.".format(name=collection_name)


class SimContainerError(Exception):
    pass


class InternalSolverNotFound(SimContainerError):
    def __init__(self, name: str) -> None:
        self.message = self.message = "Could not find iternal solver{n}".format(n=name)


class EntityTypeNotFound(SimContainerError):
    def __init__(self, name: str) -> None:
        self.message = "Could not find entity type {n}".format(n=name)
