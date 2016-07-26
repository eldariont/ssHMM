#Source: http://code.activestate.com/recipes/210459-quick-and-easy-fifo-queue-class/

class SimpleQueue:
    """A sample implementation of a First-In-First-Out
       data structure."""
    def __init__(self):
        self.in_stack = []
        self.out_stack = []

    def push(self, obj):
        self.in_stack.append(obj)

    def pop(self):
        if not self.out_stack:
            self.in_stack.reverse()
            self.out_stack = self.in_stack
            self.in_stack = []
        return self.out_stack.pop()

    def empty(self):
        return len(self.in_stack) == 0 and len(self.out_stack) == 0

    def length(self):
        return len(self.in_stack) + len(self.out_stack)

    def as_list(self):
        return self.in_stack + list(reversed(self.out_stack))