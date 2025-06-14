from typing import Type, Dict
from .file_reader import FileReader


_FILE_READER_REGISTRY: Dict[str, Type[FileReader]] = {}


def registerFileReader(alias: str):

    def decorator(cls: Type[FileReader]):

        _FILE_READER_REGISTRY[alias.upper()] = cls

        return cls

    return decorator


def fetchFileReader(alias: str, *args, **kwargs) -> FileReader:

    if (alias.upper() not in _FILE_READER_REGISTRY):

        raise ValueError()

    return _FILE_READER_REGISTRY[alias.upper()](*args, **kwargs)