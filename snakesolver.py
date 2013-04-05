#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# 
# snakesolver v0.2, 1st october 2011
#
# changelog:
#   - give all the solutions ignoring those which are equivalent by symmetry or
#     rotation
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
    
    def __init__(self, dimensions):
        """Create a volume helper.
        
        Keyword arguments:
        dimensions -- the dimensions of the target volume (e.g. [3, 3, 3])
        init_cursor -- the starting position in the volume (e.g. [0, 0, 0])
        """
        self.dimensions = dimensions
        self.path = []
        self.flags = VolumeHelper.__create_volume_flags(dimensions)
        self.init_cursor = None
    
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
    
    @staticmethod
    def __create_volume_flags(dimensions, index=0):
        """Create a multi-dimensional array filled with False values."""
        if index == len(dimensions) - 1:
            return [False] * dimensions[-1]
        return [ VolumeHelper.__create_volume_flags(dimensions, index + 1)
                 for i in xrange(dimensions[index])]
    
    def __repr__(self):
        return repr(cube_flags)

class SymmetryHelper:

    """Symmetry helper for the solver.
    
    It manages equivalence classes for equivalent points and equivalent pathes
    by symmetry and/or rotation.
    
    As symmetries and rotations only concern permutation and inversion of
    vector components (positions), then equivalences classes concern
    dimensions.
    If x and z axis are equivalent (at a specific step), we could represent the
    equivalence classes like this: [ [0, 2], [1] ]. In that case, the solver
    would only try vectors on x axis, but will ignore z axis (because it is
    equivalent).
    Then, after a move, they are not equivalent anymore, then the equivalence
    classes could be represented by [ [0], [1], [2] ].
    
    But there is a far better representation for handling quickly the
    equivalence classes: a simple array, with the same length as dimensions.
    For each dimension i, the eq_classes array contains the lower index of an
    equivalent axis:
      - i if there is no axis with lower index which is equivalent;
      - j if there is an axis j with a lower index (the lowest) which is
        equivalent.
    For example, if x and z are equivalent, then eq_classes is [ 0, 1, 0 ].
    If y and z are equivalent, then eq_classes is [ 0, 1, 1 ].
    If x and y are equivalent, then eq_classes is [ 0, 0, 2 ].
    If x, y and z are equivalent, then eq_classes is [ 0, 0, 0 ].
    If none are equivalent, then eq_classes is [ 0, 1, 2 ].
    
    With this representation, we can easily pick only 1 vector per equivalence
    class (and ignore the others): the ones with eq_classes[i] == i.
    
    eq_classes_path stores the list of eq_classes used, to restore previous
    ones when calling back(). Equivalence classes at index i are always
    computed from the equivalence classes at index i-1 (except if i == 0),
    and are always "as or more splitted".
    """

    def __init__(self, dimensions):
        self.dimensions = dimensions
        # eq_classes is always eq_classes_path[-1]
        self.eq_classes = self.__create_eq_classes_from_dimensions()
        self.eq_classes_path = [ self.eq_classes ]
    
    def __create_eq_classes_from_dimensions(self):
        """Compute the first equivalences classes from the dimensions.
        
        The dimensions which have the same length are equivalent.
        """
        eq_classes = range(len(self.dimensions))
        for i in xrange(len(self.dimensions)):
            value = self.dimensions[i]
            for j in xrange(0, i):
                # not called when i == 0
                if self.dimensions[j] == value:
                    eq_classes[i] = j
                    break
        return eq_classes
    
    def set_cursor(self, cursor):
        """Change the initial cursor position."""
        # the first item in eq_classes_path is computed from the dimensions
        # the second item is computed from the init_point
        # the next items are computed from the next moves (vectors)
        assert len(self.eq_classes_path) == 1 or\
               len(self.eq_classes_path) == 2, 'Changing cursor cannot ' +\
               'happen in the middle of a path calculation'
        if len(self.eq_classes_path) == 2:
            # this is not the first initialisation, we have to remove the
            # eq_classes associated with the previous cursor
            self.eq_classes_path.pop()
        self.cursor = cursor
        
        # make a copy to not modify the previous one (must be immutable)
        cursor_eq_classes = self.eq_classes_path[0][:]
        for i in xrange(len(cursor_eq_classes)):
            # if the equivalence class must be splitted
            # i.e. the value for the axis has changed
            if cursor_eq_classes[i] != i and not self.__eq_cmp(
                    cursor[cursor_eq_classes[i]],\
                    cursor[i], self.dimensions[i]):
                old_class = cursor_eq_classes[i]
                cursor_eq_classes[i] = i
                # If eq_classes = [ 0, 0, 0 ], and we detected that the first
                # dimension is not equivalent anymore with the two others,
                # then the new eq_classes will be [ 0, 1, 1 ]:
                # we have to check the next dimensions to change their value
                for j in xrange(i + 1, len(cursor_eq_classes)):
                    if cursor_eq_classes[j] == old_class:
                        cursor_eq_classes[j] = i
        
        self.eq_classes = cursor_eq_classes
        self.eq_classes_path.append(cursor_eq_classes)
        # At this point, eq_classes_path contains two items:
        #   - the eq_class at index 0 computed only from the dimensions;
        #   - the eq_class at index 1 computed from the initial point
        #     (which splits the equivalence classes computed from the
        #     dimensions).
        
    def get_useful_points(self, index=0, minimum=0):
        """Returns one point from each equivalence class, not more.
        
        For example, if dimensions is [3, 3, 3], then useful points are
        [[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]
        which are, respectively:
        [ a corner, an edge, a center, the core (middle of the cube) ]
        All other points are equivalent to one point in this minimal set.
        """
        if index == len(self.dimensions):
            return [[]]
        # h is "head", t is "tail"
        return ([h] + t for h in xrange(minimum,
                                        (self.dimensions[index] + 1) / 2)\
                for t in self.get_useful_points(index + 1, h))
    
    def move(self, vector):
        """Compute the new equivalent classes after a move by the vector."""
        assert self.eq_classes[vector.position] == vector.position,\
               'A move must always concern the first vector of an ' +\
               'equivalence class'
               
        # create a copy of eq_classes only if it changes, else use the same
        # instance
        has_changes = False
        position = vector.position
        new_eq_classes = self.eq_classes
        new_eq_class = None
        # the axis is now alone in its equivalence class (the vector has moved
        # the cursor on this axis, but not on the others), we have to update
        # the others (if any)
        for i in xrange(position + 1, len(self.eq_classes)):
            if self.eq_classes[i] == position:
                if not has_changes:
                    has_changes = True
                    new_eq_classes = self.eq_classes[:]
                if new_eq_class is None:
                    new_eq_class = i
                new_eq_classes[i] = new_eq_class
        self.eq_classes = new_eq_classes
        self.eq_classes_path.append(new_eq_classes)
    
    def back(self):
        """Cancel the last move."""
        self.eq_classes_path.pop()
        self.eq_classes = self.eq_classes_path[-1]
    
    def must_explore(self, i):
        """Returns True if the vector position is the first in its equivalence
        class (the others will be ignored, because equivalents)"""
        return self.eq_classes[i] == i
    
    @staticmethod
    def __eq_cmp(p1, p2, dim):
        """Comparator which tests if two point are "equivalent" on an axis.
        
        Keyword arguments:
        p1 -- projection of the point 1 on the axis
        p2 -- projection of the point 2 on the axis
        dim -- length of dimension associated to the axis
        
        p1 and p2 are equivalent if and only if:
          -    v1 == v2 (they are at the same place)
          - or v1 + v2 + 1 = dim (they are symmetrically opposite)
        """
        return p1 == p2 or p1 + p2 + 1 == dim

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
        self.symmetry_helper = SymmetryHelper(dimensions)
    
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
            for init_point in self.symmetry_helper.get_useful_points():
                # for each useful initial point (in a minimal set where no
                # points are equivalent to another)
                self.volume_helper.set_cursor(init_point)
                self.symmetry_helper.set_cursor(init_point)
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
                                     if i != previous_position and
                                        self.symmetry_helper.must_explore(i)):
                if self.volume_helper.can_move(possible_vector):
                    # if it is possible to move the cursor by the vector, then
                    # move it
                    self.volume_helper.move(possible_vector)
                    # and recursively solve
                    self.symmetry_helper.move(possible_vector)
                    for solution in self.__solve_rec(init_cursor, step + 1):
                        yield solution
                    # cancel the move to put back the state of the helper
                    self.volume_helper.back()
                    self.symmetry_helper.back()
    
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
