import pickle
from typing import Dict, List, TextIO, BinaryIO


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
    def dump(self, file_with_nb_of_objects: TextIO, pickle_file: BinaryIO):
        file_with_nb_of_objects.write(str(self.get_number_of_objects()))
        pickle.dump(self, pickle_file)

    @staticmethod
    def load(file: BinaryIO) -> "MPHF":
        return pickle.load(file)
