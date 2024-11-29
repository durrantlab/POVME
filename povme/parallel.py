from typing import Any

import multiprocessing
import platform


class MultiThreading:
    """A class for running calculations on multiple processors."""

    results: list[Any] = []

    def __init__(self, inputs, num_processors, task_class):
        """Launches a calculation on multiple processors.

        Args:
            inputs: A list, containing all the input required for the
                calculation.
            num_processors: An integer, the requested number of processors to
                use.
            task_class: An class, the class governing what calculations will be
                run on a given thread.

        Returns:
            Nothing, though the objects self.results list is populated with the
                calculation results

        """

        self.results = []

        # If it's windows, you can only use one processor.
        if num_processors != 1 and (
            platform.system().upper()[:3] == "WIN" or "NT" in platform.system().upper()
        ):
            print(
                "WARNING: Use of multiple processors is not supported in Windows. Proceeding with one processor..."
            )
            num_processors = 1

        if num_processors == 1:  # so just running on 1 processor, perhaps under windows
            single_thread = task_class()
            single_thread.total_num_tasks = len(inputs)

            single_thread.results = []
            for item in inputs:
                single_thread.value_func(item, None)

            self.results = single_thread.results

        else:  # so it actually is running on multiple processors

            cpu_count = 1
            cpu_count = multiprocessing.cpu_count()

            # first, if num_processors <= 0, determine the number of
            # processors to use programatically
            if num_processors <= 0:
                num_processors = cpu_count

            # reduce the number of processors if too many have been specified
            if len(inputs) < num_processors:
                num_processors = len(inputs)

            # if there are no inputs, there's nothing to do.
            if len(inputs) == 0:
                self.results = []
                return

            # now, divide the inputs into the appropriate number of processors
            inputs_divided = {t: [] for t in range(num_processors)}
            for t in range(0, len(inputs), num_processors):
                for t2 in range(num_processors):
                    index = t + t2
                    if index < len(inputs):
                        inputs_divided[t2].append(inputs[index])

            # now, run each division on its own processor
            running = multiprocessing.Value("i", num_processors)
            mutex = multiprocessing.Lock()

            arrays = []
            threads = []
            for _ in range(num_processors):
                athread = task_class()
                athread.total_num_tasks = len(inputs)

                threads.append(athread)
                arrays.append(multiprocessing.Array("i", [0, 1]))

            results_queue = multiprocessing.Queue()  # to keep track of the results

            processes = []
            for i in range(num_processors):
                p = multiprocessing.Process(
                    target=threads[i].runit,
                    args=(running, mutex, results_queue, inputs_divided[i]),
                )
                p.start()
                processes.append(p)

            is_running = 0  # wait for everything to finish

            while running.value > 0:
                pass
            # compile all results into one list
            for _ in threads:
                chunk = results_queue.get()
                self.results.extend(chunk)


class MultithreadingTaskGeneral:
    """A parent class of others that governs what calculations are run on each
    thread"""

    results: list[Any] = []

    def runit(self, running, mutex, results_queue, items):
        """Launches the calculations on this thread

        Args:
            running: A multiprocessing.Value object.
            mutex: A multiprocessing.Lock object.
            results_queue: A multiprocessing.Queue() object for storing the
                calculation output.
            items: A list, the input data required for the calculation.

        """

        for item in items:
            self.value_func(item, results_queue)

        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)

    def value_func(self, item, results_queue):  # so overwriting this function
        """The definition that actually does the work.

        Args:
            item: A list or tuple, the input data required for the calculation.
            results_queue: A multiprocessing.Queue() object for storing the
                calculation output.

        """

        # input1 = item[0]
        # input2 = item[1]
        # input3 = item[2]
        # input4 = item[3]
        # input5 = item[4]
        # input6 = item[5]

        # use inputs to come up with a result, some_result

        # self.results.append(some_result)

        pass


class GeneralTask:
    """A class that determines the specific calculations that will be
    performed when multi-processor support is used. Other, more
    specific classes will inherit this one."""

    results = []

    def runit(self, running, mutex, results_queue, items):
        for item in items:
            self.value_func(item, results_queue)
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)

    # this is the function that changes through inheritance
    def value_func(self, item, results_queue):
        print(item)  # here's where you do something
        self.results.append(item)  # here save the results for later compilation
