class PostProcessError(Exception):
    pass


class SubProcessTimeout(PostProcessError):

    def __init__(self, timeout: float) -> None:
        self.message = "Workers took longer than {t} seconds to finish".format(t=timeout)
        return None


class DataframeEmptyError(PostProcessError):
    def __init__(self, df_name: str ) -> None:
        self.message = "Dataframe {name} was empty".format(name=df_name)
        return None
