#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# 
# snakesolver v0.1, 27th september 2011
#
# Solver for generalized snake-cube :
# http://en.wikipedia.org/wiki/Snake_cube
# http://fr.wikipedia.org/wiki/Cube_serpent
#
# By Romain Vimont (®om)
#   rom@rom1v.com

# snake structure (list of consecutives vector norms)
# Wikipedia example
SNAKE_STRUCTURE = [ 2, 1, 1, 2, 1, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 2 ]

# size of each dimension of the target volume
VOLUME_DIMENSIONS = [ 3, 3, 3 ]

# first variable names, the next will be k1, k2, k3...
VARIABLES = [ 'x', 'y', 'z', 't' ]

class Vector:
    
    """Base vector, with only one non-zero value."""
    
    def __init__(self, position, value):
        """Create a new base vector.
        
        Examples:
          - 2 dimensions:
              Vector(0, 2) means (2, 0)
              Vector(1, 3) means (0, 3)
          - 3 dimensions:
              Vector(0, 2) means (2, 0, 0)
              Vector(1, 3) means (0, 3, 0)
          - n dimensions:
              Vector(k, i) means (0, 0, ..., 0, i at position k, 0, ..., 0, 0)
        """
        self.position = position
        self.value = value
    
    @staticmethod
    def __get_variable(position):
        """Returns the variable name associated to position.
        
        Variable names are : [ 'x', 'y', 'z', 't', 'k1', 'k2', 'k3', ... ].
        """
        if position < len(VARIABLES):
            return VARIABLES[position]
        return 'k' + str(position - len(VARIABLES) + 1)
    
    @staticmethod
    def __get_canonical(number):
        """Removes the 1 if number is 1 or -1.
        
        Used for formatting …, -3x, -2x, -x, x, 2x, 3x, …
        """
        if number == 1:
            return ''
        if number == -1:
            return '-'
        return str(number)
    
    def __repr__(self):
        return Vector.__get_canonical(self.value) +\
               Vector.__get_variable(self.position)

class VolumeHelper:
    
    """Volume helper for the solver.
    
    A volume must be considered as an "hyper-volume", it can have any number of
    dimensions. For example, a volume 3x3 is a square, while a volume 3x3x3x3
    is an hypercube.
        
    Keep a flags volume (1 boolean per "point"), indicating if the case is
    already filled.
    """
    
    def __init__(self, dimensions, init_cursor=None):
        """Create a volume helper.
        
        Keyword arguments:
        dimensions -- the dimensions of the target volume (e.g. [3, 3, 3])
        init_cursor -- the starting position in the volume (e.g. [0, 0, 0])
        """
        self.dimensions = dimensions
        self.path = []
        self.flags = VolumeHelper.__create_volume_flags(dimensions)
        self.init_cursor = None
        if init_cursor is not None:
            self.set_cursor(init_cursor)
    
    def set_cursor(self, cursor):
        """Change the initial cursor position."""
        assert len(self.path) == 0, 'Changing cursor cannot happen in the ' +\
                                    'middle of a path calculation'
        if self.init_cursor is not None:
            # set to False the flag for the old cursor
            self.__set_flag(self.init_cursor, False)
        self.init_cursor = cursor[:]
        self.cursor = self.init_cursor
        self.__set_flag(self.cursor, True)
    
    def __get_flag(self, cursor):
        """Return the flag for the specified cursor."""
        tmp = self.flags
        # dereference n times
        for i in xrange(len(cursor)):
            tmp = tmp[cursor[i]]
        # after the last iteration, tmp contains the value
        return tmp
    
    def __set_flag(self, cursor, value):
        """Change the flag for the specified cursor."""
        tmp = self.flags
        for i in xrange(len(cursor) - 1):
            tmp = tmp[cursor[i]]
        tmp[cursor[-1]] = value
    
    def can_move(self, vector):
        """Indicates if it is possible to move the cursor by the move defined
        by the vector."""
        
        # the new vector will change the cursor at position vector.position,
        # adding vector.value
        # for example, if the cursor = [1, 2, 0], and vector = Vector(2, 2),
        # then cursor will take the value [1, 2, 2]
        cursor_position_value = self.cursor[vector.position]
        new_value = cursor_position_value + vector.value
        if new_value < 0 or new_value >= self.dimensions[vector.position]:
            # outside volume
            return False
        
        # copy the cursor for not modifying the real one
        future_cursor = self.cursor[:]
        if vector.value < 0:
            sign = -1
        else:
            sign = 1
        for i in xrange(sign * vector.value):
            future_cursor[vector.position] += sign
            # if the flag at future_cursor is already True, then we cannot move
            # to this case, it is already filled
            if self.__get_flag(future_cursor):
                return False
        return True
    
    def move(self, vector):
        """Move the cursor by the vector, and updates flags and path."""
        self.path.append(vector)
        if vector.value < 0:
            sign = -1
        else:
            sign = 1
        for i in xrange(sign * vector.value):
            self.cursor[vector.position] += sign
            self.__set_flag(self.cursor, True)
    
    def back(self):
        """Cancel the last move.
        
        Used for the recursivity, for avoiding to create a new flags volume for
        each recursivity step.
        """
        vector = self.path.pop()
        if vector.value < 0:
            sign = -1
        else:
            sign = 1
        for i in xrange(sign * vector.value):
            self.__set_flag(self.cursor, False)
            self.cursor[vector.position] += -sign
    
    def all_points(self, index=0):
        """Returns all the possible points for the dimensions.
        
        For example, if dimensions is [2,2], then all_points is
        [[0,0], [0,1], [1,0], [1,1]].
        This function returns an iterator: the full list is not created.
        (replace "()" by "[]" in the return statement to construct the full
        list).
        """
        if index == len(self.dimensions):
            return [[]]
        # h is "head", t is "tail"
        return ([h] + t for h in xrange(self.dimensions[index])
                for t in self.all_points(index + 1))
    
    @staticmethod
    def __create_volume_flags(dimensions, index=0):
        """Create a multi-dimensional array filled with False values."""
        if index == len(dimensions) - 1:
            return [False] * dimensions[-1]
        return [ VolumeHelper.__create_volume_flags(dimensions, index + 1)
                 for i in xrange(dimensions[index])]
    
    def __repr__(self):
        return repr(cube_flags)
    
class SnakeCubeSolver:
    
    """Solver."""
    
    def __init__(self, dimensions, structure):
        """Create a new solver.
        
        Keyword arguments:
        dimensions -- the dimensions of the target volume
                      (for example [3, 3, 3])
        structure -- the snake structure
        """
        self.dimensions = dimensions
        self.structure = structure
        self.volume_helper = VolumeHelper(dimensions)
    
    def solve(self):
        """Solve the snake.
        
        This function returns an iterator: the full list of solutions is not
        created."""
        
        # the structure length must exactly fill the target volume
        structure_length = sum(self.structure)
        needed_length = reduce(lambda x, y: x * y, self.dimensions) - 1
        if structure_length != needed_length:
            print 'Structure has not the right length (' +\
                  str(structure_length) + ' instead of ' +\
                  str(needed_length) + ')'
        else:
            for init_point in self.volume_helper.all_points():
                # for each initial point
                self.volume_helper.set_cursor(init_point)
                # recursively solve and yield the solutions
                for solution in self.__solve_rec(init_point[:], 0):
                    yield solution
    
    def __solve_rec(self, init_cursor, step):
        """Recursively solve.
        
        Keyword arguments:
        init_cursor -- starting point, only used to put it in found solutions
        step -- recursivity depth, index of current vector in snake structure
        """
        if step == len(self.structure):
            # a new solution is found, copy the path and yield the solution
            yield init_cursor, self.volume_helper.path[:]
        else:
            if len(self.volume_helper.path) == 0:
                # empty path, use an inexistant position (-1),
                # i != previous_position will always be True
                previous_position = -1 
            else:
                previous_position = self.volume_helper.path[-1].position
            # get the vector norm for the current step
            norm = self.structure[step]
            # iterate over the next possible vectors, i.e. all vectors which
            # are orthogonal to the previous one
            for possible_vector in ( Vector(i, v)
                                     for v in [ norm, -norm ]
                                     for i in xrange(len(self.dimensions))
                                     if i != previous_position ):
                if self.volume_helper.can_move(possible_vector):
                    # if it is possible to move the cursor by the vector, then
                    # move it
                    self.volume_helper.move(possible_vector)
                    # and recursively solve
                    for solution in self.__solve_rec(init_cursor, step + 1):
                        yield solution
                    # cancel the move to put back the state of the helper
                    self.volume_helper.back()

def main():
    solver = SnakeCubeSolver(VOLUME_DIMENSIONS, SNAKE_STRUCTURE)
    # print all solutions
    for solution in solver.solve():
        print solution
    
    # print the first solutions
#    max_solutions = 5
#    solutions = solver.solve()
#    for i in xrange(max_solutions):
#        try:
#            print solutions.next()
#        except StopIteration:
#            break
    exit(0)

if __name__ == "__main__":
    main()
