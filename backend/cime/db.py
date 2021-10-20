from abc import abstractmethod


class ACimeDBO():

    def get_datasets(self):
        return self.get_datasets_by()

    @abstractmethod
    def get_datasets_by(self, **kwargs):
        pass

    @abstractmethod
    def get_dataset_by(self, **kwargs):
        pass

    @abstractmethod
    def save_dataset(self, dataset):
        pass

    @abstractmethod
    def delete_dataset_by(self, **kwargs):
        pass

    @abstractmethod
    def save_file(self, file):
        pass
