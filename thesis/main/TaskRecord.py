from time import time
from thesis.main.my_debug import message

class TaskError(Exception):
    pass

class TaskNotStartedError(TaskError):
    pass

class TaskStillRunningError(TaskError):
    pass

class TaskDoesntExistError(TaskError):
    pass

class TimeStampIsNoneError(TaskError):
    pass


class Record:

    def __init__(self, record_name):
        self.record_name = record_name
        self.child_tasks = {}

        self.running = False
        self.start_time = None
        self.end_time = None

        self.history = []
        self.info = {}

    def get_info(self):
        return self.info

    def update_child_info(self):


        for k,v in self.child_tasks.items():
            v.info.update(self.info)

            if not v.is_leaf():
                v.update_child_info()

    def start_child(self, task_name):


        if not task_name in self.child_tasks.keys():
            task = TaskRecord(task_name)
            self.child_tasks[task_name] = task

        else:
            task = self.child_tasks[task_name]

        task.info.update(self.info)

        task.start()

        return task

    def add_child(self, task):

        self.child_tasks[task.record_name] = task

    def stop_child(self, task_name):

        if task_name in self.child_tasks.keys():
            self.child_tasks[task_name].stop()

        else:
            raise TaskDoesntExistError

    def start(self):

        if self.running:
            raise TaskStillRunningError
        self.running = True

        self.start_time = time()

    def stop(self):

        if not self.running:
            message("Tried to stop task, that wasn't runnig")
            # raise TaskNotStartedError
        else:
            self.end_time = time()
            self.history.append({
                "start":self.start_time,
                "end":self.end_time,
                "info":self.get_info()
            })

            self.info = {}
            self.start_time = None
            self.end_time = None

        self.running = False

    def reset(self):
        self.stop()
        for c in self.child_tasks.values():
            c.reset()

    def get_duration(self):
        if self.running:
            raise TaskStillRunningError

        elif self.start_time is None or self.end_time is None:
            raise TimeStampIsNoneError

        else:
            return self.end_time - self.start_time

    def print(self):
        message("TIME SPENT in {n}: {t} seconds".format(t = round(self.get_duration(),2), n = self.task_name))

    def is_leaf(self):

        return len(self.child_tasks) == 0

    def gather_records(self, key_prefix = ""):

        records = {}
        for k,v in self.child_tasks.items():

            records[key_prefix + k] = v.history

            if not self.is_leaf():
                records.update(v.gather_records(key_prefix = key_prefix + k+":"))

        return records

class ClassRecord(Record):

    def get_info(self):

        return {}

class TaskRecord(Record):

    pass



