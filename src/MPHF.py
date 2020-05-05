import pickle
from typing import Dict, List


class MPHF:
    def __init__(self):
        self._object_to_id = {}
        self._id_to_object = []

    @property
    def object_to_id(self) -> Dict[object, int]:
        return self._object_to_id

    @property
    def id_to_object(self) -> List[object]:
        return self._id_to_object

    def add_object(self, object):
        if object not in self.object_to_id:
            new_id = self.get_number_of_objects()
            self.object_to_id[object] = new_id
            self.id_to_object.append(object)

    def get_id(self, object) -> int:
        return self.object_to_id[object]

    def get_object(self, object_id: int) -> object:
        return self.id_to_object[object_id]

    def get_number_of_objects(self) -> int:
        both_DS_have_the_same_length = len(self.object_to_id) == len(self.id_to_object)
        assert both_DS_have_the_same_length
        return len(self.object_to_id)

    # serialization
    # Note: not tested
    # TODO: use TextIO and BinaryIO instead?
    def save(self, file_with_nb_of_objects_filepath: str, pickle_filepath: str):
        with open(file_with_nb_of_objects_filepath, "w") as file_with_nb_of_objects_filehandler:
            file_with_nb_of_objects_filehandler.write(str(self.get_number_of_objects()))
        with open(pickle_filepath, "wb") as pickle_filehandler:
            pickle.dump(self, pickle_filehandler)

    # Note: not tested
    @staticmethod
    def load(pickle_filepath: str) -> "MPHF":
        with open(pickle_filepath, "rb") as pickle_filehandler:
            return pickle.load(pickle_filehandler)
