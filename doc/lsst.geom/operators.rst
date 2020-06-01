.. py:currentmodule:: lsst.geom

.. _lsst.geom-operators:

#############################
Operators on Point and Extent
#############################

The Point and Extent classes support many mathematical operators, but the set of available operators (and their expected behavior) isn't quite as obvious as one might think.
The table below lists all supported operators, with notes below on special cases.
For the rationale behind these operations, see :jira:`RFC-41`.

.. list-table::
   :header-rows: 1

   * - LHS
     - RHS
     - Operator
     - Result
     - Notes

   * - `PointD`
     - `PointD`
     - ``+``, ``+=``
     - Not supported
     - [#f1]_

   * - `PointD`
     - `PointD`
     - ``-``
     - `ExtentD`
     -

   * - `PointD`
     - `PointD`
     - ``-=``
     - Not supported
     - [#f2]_

   * - `PointD`
     - `PointI`
     - ``-``
     - `ExtentD`
     -

   * - `PointD`
     - `PointI`
     - ``-=``
     - Not supported
     - [#f2]_

   * - `PointD`
     - `ExtentD`
     - +, +=, -, -=
     - `PointD`
     -

   * - `PointD`
     - `ExtentI`
     - ``+``, ``+=``, ``-``, ``-=``
     - `PointD`
     -

   * - `PointI`
     - `PointD`
     - ``+``, ``+=``
     - Not supported
     - [#f1]_

   * - `PointI`
     - `PointI`
     - ``+``, ``+=``
     - Not supported
     - [#f1]_

   * - `PointI`
     - `ExtentD`
     - ``+``
     - `PointD`
     -

   * - `PointI`
     - `ExtentD`
     - ``+=``
     - Not supported
     - [#f2]_

   * - `PointI`
     - `ExtentI`
     - ``+``, ``+=``
     - `PointI`
     -

   * - `PointI`
     - `PointD`
     - ``-``
     - `ExtentD`
     -

   * - `PointI`
     - `PointD`
     - ``-=``
     - Not supported
     - [#f2]_

   * - `PointI`
     - `PointI`
     - ``-``
     - `ExtentI`
     -

   * - `PointI`
     - `PointI`
     - ``-=``
     - Not supported
     - [#f2]_

   * - `PointI`
     - `ExtentD`
     - ``-``
     - `PointD`
     -

   * - `PointI`
     - `ExtentD`
     - ``-=``
     - Not supported
     - [#f2]_

   * - `PointI`
     - `ExtentI`
     - ``-``, ``-=``
     - `PointI`
     -

   * - `ExtentD`
     - `PointD`
     - ``+``
     - `PointD`
     -

   * - `ExtentD`
     - `PointD`
     - ``+=``
     - Not supported
     - [#f2]_

   * - `ExtentD`
     - `PointD`
     - ``-``, ``-=``
     - Not supported
     - [#f1]_

   * - `ExtentD`
     - `PointI`
     - ``+``
     - `PointD`
     -

   * - `ExtentD`
     - `PointI`
     - ``+=``
     - Not supported
     - [#f2]_

   * - `ExtentD`
     - `PointI`
     - ``-``, ``-=``
     - Not supported
     - [#f1]_

   * - `ExtentD`
     - `ExtentD`
     - ``+``, ``+=``, ``-``, ``-=``
     - `ExtentD`
     -

   * - `ExtentD`
     - `ExtentI`
     - ``+``, ``+=``, ``-``, ``-=``
     - `ExtentD`
     -

   * - `ExtentI`
     - `PointD`
     - ``+``
     - `PointD`
     -

   * - `ExtentI`
     - `PointD`
     - +=
     - Not supported
     - [#f2]_

   * - `ExtentI`
     - `PointD`
     - ``-``, ``-=``
     - Not supported
     - [#f1]_

   * - `ExtentI`
     - `PointI`
     - ``+``
     - `PointI`
     -

   * - `ExtentI`
     - `PointI`
     - ``+=``
     - Not supported
     - [#f2]_

   * - `ExtentI`
     - `PointI`
     - ``-``, ``-=``
     - Not supported
     - [#f1]_

   * - `ExtentI`
     - `ExtentD`
     - ``+``, ``-``
     - `ExtentD`
     -

   * - `ExtentI`
     - `ExtentD`
     - ``+=``, ``-=``
     - Not supported
     - [#f2]_

   * - `ExtentI`
     - `ExtentI`
     - ``+``, ``-``, ``+=``, ``-=``
     - `ExtentI`
     -

   * - `ExtentD`
     - ``double``
     - ``*``, ``*=``, ``/``, ``/=``
     - `ExtentD`
     -

   * - `ExtentD`
     - ``double``
     - ``//``, ``//=``
     - Not supported
     - [#f4]_

   * - `ExtentD`
     - ``int``
     - ``*``, ``*=``, ``/``, ``/=``
     - `ExtentD`
     - [#f3]_

   * - `ExtentD`
     - ``int``
     - ``//``, ``//=``
     - Not supported
     - [#f4]_

   * - `ExtentI`
     - ``double``
     - ``*``
     - `ExtentD`
     -

   * - `ExtentI`
     - ``double``
     - ``*=``
     - Not supported
     - [#f2]_

   * - `ExtentI`
     - ``double``
     - ``/``
     - `ExtentD`
     -

   * - `ExtentI`
     - ``double``
     - ``/=``
     - Not supported
     - [#f2]_

   * - `ExtentI`
     - ``double``
     - ``//``, ``//=``
     - Not supported
     - [#f4]_

   * - `ExtentI`
     - ``int``
     - ``*``, ``*=``
     - `ExtentI`
     -

   * - `ExtentI`
     - ``int``
     - ``/``
     - `ExtentD`
     -

   * - `ExtentI`
     - ``int``
     - ``/=``
     - Not supported (Python), `ExtentI` (C++)
     - [#f2]_

   * - `ExtentI`
     - ``int``
     - ``//``, ``//=``
     - `ExtentI`
     - [#f5]_

   * - ``double``
     - `ExtentD`
     - ``*``
     - `ExtentD`
     -

   * - ``double``
     - `ExtentI`
     - ``*``
     - `ExtentD`
     -

   * - ``int``
     - `ExtentD`
     - ``*``
     - `ExtentD`
     - [#f3]_

   * - ``int``
     - `ExtentI`
     - ``*``
     - `ExtentI`
     -

.. [#f1] Operation is not geometrically meaningful.
.. [#f2] This is an in-place operator that would require the LHS type to change.
         That would actually be possible to implement in Python, but its behavior would be confusing.
.. [#f3] This operator is not implemented directly in either C++ or Python, but is largely supported by the fact that an overload that takes ``double`` will also accept ``int`` (but may yield different answers for extremely large integers that cannot be represented exactly as ``double``\s).
.. [#f4] The ``//`` operator applies only to integer types.
.. [#f5] This Python-only operation does not always produce the same result as regular division of integers in C++, because Python specifies that ``a//b`` is equivalent to ``floor(a/b)``, while C++ specifies that it should be equivalent to ``int(a/b)``.
         Note that ``floor`` rounds negative numbers down and ``int`` rounds them up.
