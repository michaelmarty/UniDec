"""
Implements 'queue' data structure in Python
"""
import os


class Node:
    def __init__(self, data):
        self.data = data
        self.next = None

    def __str__(self):
        return self.data


class Queue:
    def __init__(self):
        self.__head = None
        self.__tail = None

    def __len__(self):
        if self.is_empty:
            return 0
        curr = self.__head
        counter = 0
        while curr:
            counter += 1
            curr = curr.next
        return counter

    def __iter__(self):
        curr = self.__head
        while curr:
            yield curr
            curr = curr.next

    def __contains__(self, item):
        curr = self.__head
        while curr:
            if curr.data == item:
                return True
            curr = curr.next
        return False

    def __str__(self):
        delim = '->'
        my_str = ''
        if self.is_empty:
            return 'None'
        curr = self.__head
        while curr is not None:
            my_str += f'{curr}{delim}'
            curr = curr.next
        if curr is None:
            my_str += f'None'
        return my_str

    def __set_tail_none(self):
        self.__tail = None

    @property
    def is_empty(self):
        if self.__head is None and self.__tail is None:
            return True
        else:
            return False

    def enqueue(self, val):
        new_node = Node(val)
        if self.is_empty:
            self.__head = self.__tail = new_node
        else:
            self.__tail.next = new_node
            self.__tail = new_node

    def dequeue(self):
        if self.is_empty:
            return None
        dequeued_node = self.__head
        self.__head = self.__head.next
        if self.__head is None:
            self.__tail = None
        return dequeued_node

    def reverse(self):

        pass

    def peek(self):
        return self.__head

    def close(self):
        return self.__tail


def main():
    print('Checking Methods')
    print()
    my_queue = Queue()
    print(f'Checking __str__() methods:\nEmpty Handler:\n{my_queue}')
    my_queue.enqueue('First Person')
    my_queue.enqueue('Second Person')
    my_queue.enqueue('Third Person')
    my_queue.enqueue('Fourth Person')
    print(f'Full Handler:\n{my_queue}')
    print()
    print('Testing __iter__() methods')
    for node in my_queue:
        print(node)
    print(my_queue)
    print()
    print('Testing Property Methods')
    print(my_queue.peek())
    print(my_queue.close())
    print(my_queue.is_empty)
    print()
    print('Testing dequeue:')
    print(f'Removing: {my_queue.peek()}')
    my_queue.dequeue()
    for i in range(len(my_queue)):
        print(f'Removing: {my_queue.dequeue()}')
    print(my_queue)
    print('Making sure __set...() method is not reachable:')
    try:
        my_queue.__set_tail_none()
        print("It doesn't work")
    except AttributeError:
        print('It works')
    job_list = [print('my mom was here'), "was here the sequence"]
    for job in job_list:
        my_queue.enqueue(job)
    while not my_queue.is_empty:
        my_queue.dequeue()
    print()
    any(my_queue.enqueue(file) for file in os.listdir('Yeast_Genome'))
    print(f'{my_queue.peek()}')
    try:
        os.listdir('Fake Folder')
    except FileNotFoundError:
        print('file not found')


if __name__ == '__main__':
    main()
