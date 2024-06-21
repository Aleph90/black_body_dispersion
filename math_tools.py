"""Supporting math tools for ray tracing.

This module implements some useful classes for ray tracing, including:
    - Vect3D: 3-dimensional vectors featuring linear algebra and Euclidean geometry.
    - Matrix3by3: 3 by 3 matrices with focus on matrix-vector multiplication.
    - Ray: Oriented half-lines in 3-dimensional space, determined by an origin and direction.
    - RGBColor: Subclass of Vect3D to describe color.
    - SumFun: A class to represent math functions obtained as sum of elementary components.
"""

from math import sqrt


# Geometry and linear algebra
class Vect3D:
    """A custom-made class for 3-dimensional linear objects.

    The class features access to components as ``x``, ``y``, and ``z`` as well as by index, and all the basic
    operations from linear algebra and Euclidean geometry.
    """

    def __init__(self, x: float, y: float, z: float):
        """ Initializer for the Vect3D class.

        The three components are organized into a python list.

        :param float x: ``x``-component of the new vector.
        :param float y: ``y``-component of the new vector.
        :param float z: ``z``-component of the new vector.
        """

        self._entries = [x, y, z]

    @classmethod
    def cast(cls, v: 'Vect3D') -> 'Vect3D':
        """Alternative initializer.

        Allows to initialize a vector from another vector. This method will be used later
        by subclasses to ensure that algebraic operations stay within the subclass.

        :param Vect3D v: ``Vect3D`` to be cast.
        :return: A new vector with the same entries as ``v``, instance of the same (sub-)class
            as the object to which the method is applied.
        """

        return cls(*v.entries)

    @classmethod
    def null(cls):
        """The null vector, implemented as the return of a class method.

        :return: A new ``Vect3D`` object with all components equal to ``0``, instance of the same
            (sub-)class as the object to which the method is applied.
        """

        return cls(0, 0, 0)

    def __eq__(self, other: 'Vect3D') -> bool:
        """Equality operator for ``Vect3D`` objects.

        :param Vect3D other: ``Vect3D`` object to compare ``self`` with.
        :return: ``True`` if the entries of the two vectors are equal as lists, ``False`` otherwise.
        """

        return self._entries == other.entries

    def __copy__(self) -> 'Vect3D':
        """Copy operator.

        Use the cast method to make a copy of ``self``. This implementation ensures that
        the method returns an object of the same class as ``self`` even when applied to
        an instance of a subclass.

        :return: A new vector with the same class and components as ``self``.
        """

        return self.cast(self)

    def copy(self) -> 'Vect3D':
        """Explicit copy method for the ``Vect3D`` class.

        :return: A new vector with the same class and components as ``self``.
        """
        return self.__copy__()

    # Representations and type conversion
    def __repr__(self) -> str:
        """Representation operator for the ``Vect3D`` class.

        :return: A string of the form ``'Vect3D(x, y, z)'``, with ``x``, ``y``, and ``z`` the components of ``self``.
        """

        return f'Vect3D{tuple(self._entries)}'

    def __str__(self) -> str:
        """Conversion operator ino ``str`` for the ``Vect3D`` class.

        :return: A string of the form ``(x, y, z)``, where ``x``, ``y``, and ``z`` are the components of ``self``.
        """

        return f'{tuple(self._entries)}'

    def __bool__(self) -> bool:
        """Boolean operator for the ``Vect3D`` class.

        :return: ``True`` if at least one component of ``self`` is non-zero, ``False`` otherwise.
        """
        return bool([e for e in self._entries if e])

    # Access
    @property
    def entries(self) -> list[float]:
        """Entries of a ``Vect3D`` object (read only).

        :return: The list of entries of ``self``.
        """

        return self._entries

    @property
    def x(self) -> float:
        """``x``-component.

        :return: The value of the ``x``-component of ``self``.
        """

        return self._entries[0]

    @property
    def y(self) -> float:
        """``y``-component.

        :return: The value of the ``y``-component of ``self``.
        """

        return self._entries[1]

    @property
    def z(self) -> float:
        """``z``-component.

        :return: The value of the ``z``-component of ``self``.
        """

        return self._entries[2]

    def __getitem__(self, key: int) -> float:
        """Getter by index.

        :param int key: An index between 0 and 2 (inclusive).
        :return: The value of the ``key``-th component of ``self``.
        """

        return self._entries[key]

    @x.setter
    def x(self, value: float) -> None:
        """``x``-component (setter).

        :param float value: Value to be assigned to the ``x``-component of ``self``.
        :return: ``None``.
        """

        self._entries[0] = value

    @y.setter
    def y(self, value: float) -> None:
        """``y``-component (setter).

        :param float value: Value to be assigned to the ``y``-component of ``self``.
        :return: ``None``.
        """

        self._entries[1] = value

    @z.setter
    def z(self, value: float) -> None:
        """``z``-component (setter).

        :param float value: Value to be assigned to the ``z``-component of ``self``.
        :return: ``None``.
        """

        self._entries[2] = value

    def __setitem__(self, key: int, value: float) -> None:
        """Setter by index.

        :param int key: Index of the entry to be assigned, between 0 and 2 (inclusive).
        :param float value: Value to be assigned to the ``key``-th entry.
        :return: ``None``.
        """

        self._entries[key] = value

    # Linear algebra
    def __add__(self, other: 'Vect3D') -> 'Vect3D':
        """Component-wise addition for ``Vect3D`` objects.

        :param Vect3D other: The second summand.
        :return: A new ``Vect3D`` object representing the component-wise sum of ``self`` and ``other``.
        """
        return Vect3D(*[self._entries[i] + other[i] for i in range(3)])

    def __sub__(self, other: 'Vect3D') -> 'Vect3D':
        """Component-wise subtraction for ``Vect3D`` objects.

        :param Vect3D other: The ``Vect3D`` object to be subtracted from ``self``.
        :return: A new ``Vect3D`` object representing the component-wise difference of ``self`` and ``other``.
        """

        return self + (-other)

    def __rmul__(self, scalar: float) -> 'Vect3D':
        """Component-wise scalar-by-vector multiplication.

        :param float scalar: The scaling factor.
        :return: A new ``Vect3D`` object representing ``self`` rescaled by the specified factor.
        """

        return Vect3D(*[scalar * entry for entry in self._entries])

    def __mul__(self, scalar: float) -> 'Vect3D':
        """Component-wise vector-by-scalar multiplication.

        Based on scalar-by-vector implementation.

        :param float scalar: The scaling factor.
        :return: A new ``Vect3D`` object representing ``self`` rescaled by the specified factor.
        """

        return scalar * self

    def __truediv__(self, scalar: float) -> 'Vect3D':
        """Component-wise division of a vector by a scalar.

        Based on scalar-by-vector multiplication.

        :param float scalar: The number by which to divide ``self``.
        :return: A new ``Vect3D`` object representing ``self`` rescaled by the *inverse* of the specified factor.
        """

        return (1/scalar) * self

    def __neg__(self) -> 'Vect3D':
        """Component-wise negative for ``Vect3D`` objects.

        :return: A new ``Vect3D`` object whose components are the opposite of those of ``self``.
        """

        return (-1) * self

    # Euclidean structure
    def dot_with(self, other: 'Vect3D') -> float:
        """Dot product for ``Vect3D`` objects.

        :param Vect3D other: The second factor in the dot product.
        :return: The dot product of ``self`` and ``other``.
        """

        return sum([self._entries[i] * other[i] for i in range(3)])

    def norm_squared(self) -> float:
        """Norm squared of a vector.

        Based on ``dot_with``.

        :return: Norm squared of ``self``.
        """

        return self.dot_with(self)

    def norm(self) -> float:
        """Norm of a vector.

        Based on ``norm_squared``.

        :return: Norm of ``self``.
        """

        return sqrt(self.norm_squared())

    def normalized(self) -> 'Vect3D':
        """Unit vector pointing in the same direction as ``self``.

        The returned vector is obtained by rescaling ``self`` by the inverse of its norm.
        No effect on ``self``.

        :return: A new ``Vect3D`` object representing a unit vector pointing in the same direction as ``self``.
        """

        return (1/self.norm()) * self

    def cross_with(self, other: 'Vect3D') -> 'Vect3D':
        """Cross product for the ``Vect3D`` class.

        :param Vect3D other: Second factor of the cross product.
        :return: A new ``Vect3D`` object representing the cross product of ``self`` and ``other``.
        """

        return Vect3D(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)

    def reflect(self, normal: 'Vect3D'):
        """Reflection about a plane of given normal.

        :param ``Vect3D`` normal: Vector normal to the plane of reflection.
        :return: A new vector representing the reflection of ``self`` with respect to the plane of given normal.
        """

        return self - 2 * self.dot_with(normal)/normal.norm_squared() * normal

    def refract(self, normal: 'Vect3D', eta_in: float, eta_out: float):
        """Refraction when crossing a plane of given normal.

        According to Snell's law, if a beam of light crosses a surface between a medium of refractive index
        ``eta_in`` to one of index ``eta_out``, then the direction of propagation changes according to the relation
            `eta_in * sin_in = eta_out * sin_out`
        where ``sin_in`` and ``sin_out`` express the sine of the incoming and outgoing directions of propagation
        relative to the normal direction of the surface, respectively.
        Given the angle of the incoming beam, refraction may not always be possible, as when
        ``eta_in * sin_in / eta_out > 1`` there is no possible angle whose sine satisfies the relation above. In
        that case, the beam is *reflected* instead.

        :param Vect3D normal: Vector normal to the plane.
        :param eta_in: Refractive index of the current material.
        :param eta_out: Refractive index of the material beyond the plane.
        :return: A new (unit) vector pointing in the direction after refraction if possible,
            or after reflection otherwise.
        """

        r_in = self.normalized()
        n = normal.normalized()
        d = eta_out**2 - eta_in**2 * (1 - dot(r_in, n))**2
        if d >= 0:
            # Refraction is possible! The refracted vector is a unit vector which is
            # a linear combination of r_in and n. The exact coefficients are obtained
            # by setting up and solving appropriate equations.
            alpha = - eta_in * dot(r_in, n)/sqrt(d)
            r_out = normalized((1-alpha)*dot(r_in, n)*n + alpha*r_in)
        else:
            r_out = r_in.reflect(n)
        return r_out


class Matrix3by3:
    """Custom class for 3 by 3 matrices.

    The primary use of this class will be for matrix-vector multiplication, implemented
    as a linear combination of columns. As such, a matrix is implemented as a triple of
    `Vect3D` objects representing its _columns_, contrary to the common use of
    representing a matrix as a list of rows.
    """

    def __init__(self, c1: Vect3D, c2: Vect3D, c3: Vect3D):
        """Initializer of the `Matrix3by3 class.

        :param Vect3D c1: First column.
        :param Vect3D c2: Second column.
        :param Vect3D c3: Third column.
        """

        self._cols = [c1.copy(), c2.copy(), c3.copy()]

    def __eq__(self, other: 'Matrix3by3') -> bool:
        """Equality operator for ``Matrix3by3`` objects.

        :param Matrix3by3 other: ``Matrix3by3`` object to be compared with ``self``.
        :return: ``True`` if all three columns of both matrices are equal, ``False`` otherwise.
        """
        return self._cols == other._cols

    def __copy__(self) -> 'Matrix3by3':
        """Shallow copy operator for ``Matrix3by3`` objects.

        :return: A new matrix with the same columns as ``self``.
        """

        return Matrix3by3(*[c for c in self._cols])

    def __deepcopy__(self) -> 'Matrix3by3':
        """Deep copy operator for ``Matrix3by3`` objects.

        :return: A new matrix whose columns are copies of those of ``self``.
        """

        return Matrix3by3(*[c.copy() for c in self._cols])

    # Representation and type conversion
    def __repr__(self) -> str:
        """Representation operator for the ``Matrix3by3`` class.

        :return: A string representing ``self``.
        """

        return f'Matrix3by3({", ".join([c.__repr__() for c in self._cols])})'

    def __str__(self) -> str:
        """String conversion operator for the ``Matrix3by3`` class.

        :return: A string displaying ``self`` as a list of rows.
        """

        return f'[{[c.x for c in self._cols]},\n{[c.y for c in self._cols]},\n{[c.z for c in self._cols]}]'

    def __bool__(self) -> bool:
        """Boolean operator for the ``Matrix3by3`` class.

        :return: ``True`` if ``self`` contains any one non-zero entry, ``False`` otherwise.
        """

        return bool([c for c in self._cols if c])

    # Access
    def __getitem__(self, pos: int | tuple[int, int]) -> Vect3D | float:
        """Getter operator for an index or pair of indices.

        The expected order of the indices is as row, then column.
        Since the data of ``self`` is implemented as a list of columns,
        the order of the indices needs to be reversed for access.

        :param int | tuple[int, int] pos: Index or pair of indices of the position to access.
        :return: Column in position ``pos`` if it is an integer, or value of the entry in position ``pos``.
        """
        if isinstance(pos, int):
            return self._cols[pos]
        i, j = pos
        return self._cols[j][i]

    def __setitem__(self, pos: tuple[int, int], val: float) -> None:
        """Setter operator by pair of indices.

        The expected order of the indices is as row, then column.
        Since the data of ``self`` is implemented as a list of columns,
        the order of the indices needs to be reversed for access.

        :param tuple[int, int] pos: Pair of indices of the position to access.
        :param float val: Value to be assigned to the entry in the specified position.
        :return: ``None``.
        """

        i, j = pos
        self._cols[j][i] = val

    # Algebra
    def __rmul__(self, scalar: float) -> 'Matrix3by3':
        """Multiplication by a scalar from the left.

        :param float scalar: The factor by which to rescale ``self``.
        :return: A new ``Matrix3by3`` object representing the rescaled matrix.
        """

        return Matrix3by3(*[scalar * col for col in self._cols])

    def __mul__(self, other: 'float | Vect3D | Matrix3by3'):
        """Multiplication operator for the ``Matrix3by3`` class.

        The right factor may be:
            - A scalar (i.e. float): Use ``__rmul__()``;
            - A vector (i.e. ``Vect3D`` object): Implement product as a linear combination of columns;
            - A matrix (i.e. ``Matrix3by3`` object): Form a new matrix from the matrix-by-vector products.

        If all the columns of the matrix are instances of some subclass of
        `Vect3D`, then the result will also be an instance of that class.

        :param Union[float, Vect3D, Matrix3by3] other: Right factor.
        :return: A new ``Matrix3by3`` or ``Vect3D`` representing the product.
        """

        if isinstance(other, float):
            product = other * self
        elif isinstance(other, Vect3D):
            product = sum([other[i]*self._cols[i] for i in range(3)], start=self._cols[0].null())
        elif isinstance(other, Matrix3by3):
            product = Matrix3by3(*[self * col for col in other._cols])
        else:
            raise TypeError('What are you trying to multiply this matrix by?')
        return product

    def det(self) -> float:
        """Determinant function for the ``Matrix3by3`` class.

        :return: The determinant of ``self``.
        """
        return dot(self._cols[0], cross(self._cols[1], self._cols[2]))
    
    def gram_schmidt(self, tolerance: int = 3):
        """Orthonormalize the frame formed by the columns of ``self`` using the Gram-Schmidt process.

        The Gram-Schmidt process produces an orthonormal basis out of a given frame.
        The idea is to remove from each vector of the frame its projection onto the
        span of the previous elements in the basis, and then normalizing it. Once
        the first ``i`` vectors are processed, finding the orthogonal projection of
        the next element onto their span boils down to a sequence of dot products.

        If the given vectors are not linearly independent then one of them will lie
        in the span of the previous elements. In that case, that vector reduces to
        zero after removing the orthogonal projection, and the normalization will
        raise an error. Since detecting linear dependence is so simple, this
        implementation does not check for linear independence before starting the
        process. Instead, it checks that the norm of each vector is larger than a
        certain cutoff just before normalizing it and raises an error otherwise.

        :param int tolerance: Number of digits of precision when checking for vanishing of norms.
        :return: A new ``Matrix3by3`` formed by the vectors in the orthogonalized frame.
        """

        ortho_vectors = []
        for col in self._cols:
            v = col.copy()
            for e in ortho_vectors:
                v -= dot(v, e) * e
            n = v.norm()
            if n <= 10**(-tolerance):
                raise ValueError('Gram-Schmidt fail: Columns not linearly independent.')
            v /= n
            ortho_vectors.append(v)
        return Matrix3by3(*ortho_vectors)


def dot(v: Vect3D, w: Vect3D) -> float:
    """Dot product of two vectors as a stand-alone function.

    :param Vect3D v: First "factor".
    :param Vect3D w: Second "factor".
    :return: The dot product of ``v`` and ``w``.
    """

    return v.dot_with(w)


def norm_squared(v: Vect3D) -> float:
    """Norm squared of a vector as a stand-alone function.

    :param Vect3D v: The vector.
    :return: The norm squared of ``v``.
    """

    return v.norm_squared()


def norm(v: Vect3D) -> float:
    """The norm of a vector as a stand-alone function.

    :param Vect3D v: The vector.
    :return: The norm of ``v``.
    """

    return v.norm()


def normalized(v: Vect3D) -> Vect3D:
    """Normalization of a vector as a stand-alone function.

    :param Vect3D v: The vector.
    :return: The normalization of ``v``.
    """

    return v.normalized()


def cross(v: Vect3D, w: Vect3D) -> Vect3D:
    """Cross product as a stand-alone function.

    :param Vector3D v: Left factor.
    :param Vector3D w: Right factor.
    :return: The cross product of ``v`` and ``w``.
    """

    return v.cross_with(w)


def det(mat: Matrix3by3) -> float:
    """Determinant of a 3 by 3 matrix as a stand-alone function.

    :param Matrix3by3 mat: The matrix.
    :return: The determinant of ``mat``.
    """
    return mat.det()


def is_basis(u: Vect3D, v: Vect3D, w: Vect3D, tolerance: int = 6):
    """Check if a given triple of ``Vect3D`` objects is a basis.

    Linear independence is tested via the determinant. Since the latter may
    be subject to floating point errors, the vanishing is checked up to an
    adjustable number of significant digits.

    :param Vect3D u: First vector.
    :param Vect3D v: Second vector.
    :param Vect3D w: Third vector.
    :param int tolerance: Number of digits of precision when checking for vanishing of the determinant.
    :return: ``True`` if the given vectors are linearly independent (up to the specified error), ``False`` otherwise.
    """

    return abs(det(Matrix3by3(u, v, w))) > 10**(-tolerance)


def is_positive_frame(u: Vect3D, v: Vect3D, w: Vect3D, tolerance: int = 6):
    """Check if a given triple of ``Vect3D`` objects is a positively oriented basis.

    This is done by checking that the determinant of the corresponding matrix
    is positive, or rather greater than an adjustable cutoff to account for
    floating point errors.

    :param Vect3D u: First vector.
    :param Vect3D v: Second vector.
    :param Vect3D w: Third vector.
    :param int tolerance: Digits of precision to be considered when checking for vanishing of the determinant.
    :return: ``True`` if the three vectors form a positively oriented frame, ``False`` otherwise.
    """

    return det(Matrix3by3(u, v, w)) > 10**(-tolerance)


class Ray:
    """A class for rays in 3 dimensions.

    A ray represents an oriented (half) line emanating from a given point,
    its *origin*, along a given vector, its *direction*. In ray tracing,
    this kind of object is used to model the path of light emitted by a
    source or collected by a receptor.

    Given the data of the origin and direction of a ray, it is straightforward
    to obtain a parameterization of the corresponding (half) line and
    therefore reduce the problem of finding intercept with objects to
    a functional equation in one real parameter.

    Attributes
    ----------
    origin : Vect3D
        The coordinate vector of the point from which the ray emanates.
    direction : Vect3D
        The vector along which the ray emanates from its origin.
    """

    def __init__(self, origin: Vect3D, direction: Vect3D):
        """Initializer of the ``Ray`` class.

        A ray is constructed out of the data of an origin and direction,
        which are copied and stored as instance attributes.

        :param Vect3D origin: Base point of the ray.
        :param Vect3D direction: Direction vector of the ray.
        """

        self.origin = origin.copy()
        self.direction = direction.copy()

    def __eq__(self, other: 'Ray') -> bool:
        """Equality operator for the ``Ray`` class.

        :param Ray other: The ray with which to compare ``self``.
        :return: ``True`` if ``self`` and ``other`` have equal origins and directions, ``False`` otherwise.
        """

        return self.origin == other.origin and self.direction == other.direction

    def __copy__(self) -> 'Ray':
        """Copy operator for the ``Ray`` class.

        :return: A new ``Ray`` object with the same origin and direction as ``self``.
        """

        return Ray(self.origin, self.direction)

    def copy(self) -> 'Ray':
        """Explicit copy operator for the ``Ray`` class.

        :return: A new ``Ray`` object with the same origin and direction as ``self``.
        """

        return self.__copy__()

    def __repr__(self) -> str:
        """Representation operator for the ``Ray`` class.

        :return: A string representation of ``self``.
        """

        return f'Ray({str(self.origin)}, {str(self.direction)})'

    def __str__(self) -> str:
        """String conversion for the ``Ray`` class.

        :return: A string expressing ``self`` in the form of an ordered pair
            consisting of origin and direction.
        """
        
        return f'({str(self.origin)}, {str(self.direction)})'

    def point_at(self, t: float | None) -> Vect3D | None:
        """Determine the point on the parameterized line at time ``t``.
        
        :param float | None t: Time, optional.
        :return: The coordinate vector of the point on the parameterized
            line corresponding to the value ``t`` of the parameter if it is
            specified, ``None`` otherwise.
        """
        
        return self.origin + t * self.direction if t is not None else None


class SumFun:
    """A class of real-valued functions (in the math sense) that come as a sum of more fundamental ones.

    Some of the important functions we will use later, such as the color matching functions
    and Sellmeier's law, occur as sums of fundamental constituents. It is implicitly assumed
    that all these summands should all have consistent input and return ``float`` type objects.
    """

    def __init__(self, *components):
        """Initializer for the SumFun class.

        The components to be summed are stored as a list in the protected attribute ``_components``.
        The resulting object represents the mathematical function obtained as the sum of all the
        functions in the list.

        No checks are carried out on the parameters ``components`` at this stage. However, it is
        implicitly assumed that they should all have input types consistent with each other and
        return ``float`` type objects.

        :param components: The functions of which the new instance represents the sum.
        """

        self._components = list(components)

    # Representation and type conversion
    def __repr__(self):
        """Representation operator for the ``SumFun`` class.

        :return: A string representation of the form ``'SumFun(f_1, f_2, ..., f_n)'``.
        """

        return f'SumFun{tuple(self._components)}'

    def __str__(self):
        """String conversion operator for the ``SumFun`` class.

        :return: A string representation of the form ``'f_1 + f_1 + ... + f_n'``.
        """

        return ' + '.join([str(comp) for comp in self._components])

    # Function call
    def __call__(self, *args) -> float:
        """Call operator for the ``SumFun`` class.

        Instances of the ``SumFun`` class represent functions and should
        therefore be callable. The collection of parameters ``args`` is
        assumed to be of the common input type of all the summands in
        `self._components`. The operator simply evaluates all component
        functions at the given ``args`` and returns the sum of the outputs.

        :param args: Arguments to be passed to all components of ``self``.
        :return: The sum of the values of all components at ``args``.
        """

        return sum([component(*args) for component in self._components])

    # Algebra
    def __add__(self, other: 'SumFun') -> 'SumFun':
        """Addition operator for the ``SumFun`` class.

        Since instances of ``SumFun`` are essentially encoded as lists
        of component functions, their sum is implemented as a new
        object formed by the concatenated lists of components.

        :param SumFun other: The ``SumFun`` function to be added to ``self``.
        :return: A new function representing the sum of ``self`` and ``other``.
        """

        return SumFun(*self._components + other._components)
