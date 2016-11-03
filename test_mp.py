#!/usr/bin/env python3
# --*-- coding:utf-8 --*--
from random import randint
from time import ctime
from time import sleep
import queue
import threading


class MyTask(object):
    def __init__(self, name):
        self.name = name
        self._work_time = randint(1, 5)

    def work(self):
        print("Task %s is start : %s, sleep time= %d" % (self.name, ctime(), self._work_time))
        sleep(self._work_time)
        print("Task %s is end : %s" % (self.name, ctime()))


class MyThread(threading.Thread):
    def __init__(self, my_queue):
        self.my_queue = my_queue
        super(MyThread, self).__init__()

    def run(self):participle
        while True:
            if self.my_queue.qsize() > 0:
                self.my_queue.get().work()
            else:
                break


def print_split_line(num=30):
    print("*" * num)

if __name__ == "__main__":
    print_split_line()
    queue_length = 6
    my__queue = queue.LifoQueue(queue_length * 3)
    threads = []

    for i in range(queue_length * 3):
        mt = MyTask("tk_" + str(i))
        my__queue.put_nowait(mt)

    for i in range(queue_length):
        mtd = MyThread(my__queue)
        threads.append(mtd)

    for i in range(queue_length):
        threads[i].start()

    for i in range(queue_length):
        threads[i].join()

    print_split_line()
