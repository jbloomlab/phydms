class Base:

    def __init__(self, s):
        print("in Base __init__: {0}".format(s))
        self.x(s)

    def x(self, s):
        print("in Base x: {0}".format(x))


class Child(Base):
    
    def __init__(self, s):
        print("in Child __init__: {0}".format(s))
        super(Child, self).__init__(s=s)
        self.x(s)

    def x(self, s):
        print("in Child x: {0}".format(s))


c = Child('Lauren')


