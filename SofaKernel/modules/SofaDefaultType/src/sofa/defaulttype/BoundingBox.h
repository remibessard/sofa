/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_DEFAULTTYPE_BOUNDINGBOX_H
#define SOFA_DEFAULTTYPE_BOUNDINGBOX_H

#include <utility> /// for std::pair
#include <sofa/defaulttype/config.h>
#include <sofa/defaulttype/Vec.h>

namespace sofa
{
namespace defaulttype
{

class SOFA_DEFAULTTYPE_API BoundingBox
{

public:
    typedef std::pair< Vector3, Vector3 > bbox_t;

    BoundingBox();
    /// Define using the endpoints of the main diagonal
    BoundingBox(const Vector3& minBBox, const Vector3& maxBBox);
    BoundingBox(const bbox_t& bbox);
    /// Define using xmin, xmax, ymin, ymax, zmin, zmax in this order
    BoundingBox(SReal xmin, SReal xmax, SReal ymin, SReal ymax, SReal zmin, SReal zmax );
    /// Define using xmin, xmax, ymin, ymax, zmin, zmax in this order
    BoundingBox(const Vec6f& bbox);
    /// Define using xmin, xmax, ymin, ymax, zmin, zmax in this order
    BoundingBox(const Vec6d& bbox);

    static BoundingBox neutral_bbox();

    operator bbox_t() const;

    void invalidate();
    bool isValid() const;
    bool isFlat()  const;
    bool isNegligeable() const; // !valid || flat
    bool isNull()  const;

    SReal* minBBoxPtr();
    SReal* maxBBoxPtr();
    const SReal* minBBoxPtr() const;
    const SReal* maxBBoxPtr() const;
    const Vector3&  minBBox() const;
    const Vector3&  maxBBox() const;
    Vector3& minBBox();
    Vector3& maxBBox();

    bool contains( const sofa::defaulttype::Vector3& point) const;
    bool contains( const BoundingBox& other) const;

    bool intersect( const BoundingBox& other) const;
    void intersection( const BoundingBox& other);

    void include( const sofa::defaulttype::Vector3& point);
    void include( const BoundingBox& other);

    void inflate( SReal amount );

    BoundingBox getIntersection( const BoundingBox& other ) const;
    BoundingBox getInclude( const sofa::defaulttype::Vector3& point ) const;
    BoundingBox getInclude( const BoundingBox& other ) const;
    BoundingBox getInflate( SReal amount ) const;

    friend std::ostream& operator << ( std::ostream& out, const BoundingBox& bbox)
    {
        out << bbox.minBBox() << " " <<  bbox.maxBBox();
        return out;
    }

    friend std::istream& operator >> ( std::istream& in, BoundingBox& bbox)
    {
        in >> bbox.minBBox() >> bbox.maxBBox();
        return in;
    }


protected:
    bbox_t bbox;
};


template <typename TReal>
class TBoundingBox : public BoundingBox
{
public:
    TBoundingBox(const TReal* minBBoxPtr, const TReal* maxBBoxPtr)
        :BoundingBox(Vector3(minBBoxPtr),Vector3(maxBBoxPtr))
    {
    }

};



/// bounding rectangle
class SOFA_DEFAULTTYPE_API BoundingBox2D
{

public:
    typedef std::pair< Vector2, Vector2 > bbox_t;

    BoundingBox2D();
    /// Define using the endpoints of the main diagonal
    BoundingBox2D(const Vector2& minBBox, const Vector2& maxBBox);
    BoundingBox2D(const bbox_t& bbox);
    /// Define using xmin, xmax, ymin, ymax in this order
    BoundingBox2D(SReal xmin, SReal xmax, SReal ymin, SReal ymax );
    /// Define using xmin, xmax, ymin, ymax in this order
    BoundingBox2D(const Vec4f& bbox);
    /// Define using xmin, xmax, ymin, ymax in this order
    BoundingBox2D(const Vec4d& bbox);

    static BoundingBox2D neutral_bbox();

    operator bbox_t() const;

    void invalidate();
    bool isValid() const;
    bool isFlat()  const;
    bool isNegligeable() const; // !valid || flat
    bool isNull()  const;

    SReal* minBBoxPtr();
    SReal* maxBBoxPtr();
    const SReal* minBBoxPtr() const;
    const SReal* maxBBoxPtr() const;
    const Vector2&  minBBox() const;
    const Vector2&  maxBBox() const;
    Vector2& minBBox();
    Vector2& maxBBox();

    bool contains( const sofa::defaulttype::Vector2& point) const;
    bool contains( const BoundingBox2D& other) const;

    bool intersect( const BoundingBox2D& other) const;
    void intersection( const BoundingBox2D& other);

    void include( const sofa::defaulttype::Vector2& point);
    void include( const BoundingBox2D& other);

    void inflate( SReal amount);

    BoundingBox2D getIntersection( const BoundingBox2D& other ) const;
    BoundingBox2D getInclude( const sofa::defaulttype::Vector2& point ) const;
    BoundingBox2D getInclude( const BoundingBox2D& other ) const;
    BoundingBox2D getInflate( SReal amount ) const;

    friend std::ostream& operator << ( std::ostream& out, const BoundingBox2D& bbox)
    {
        out << bbox.minBBox() << " " <<  bbox.maxBBox();
        return out;
    }

    friend std::istream& operator >> ( std::istream& in, BoundingBox2D& bbox)
    {
        in >> bbox.minBBox() >> bbox.maxBBox();
        return in;
    }


protected:
    bbox_t bbox;
};



/// bounding interval
class SOFA_DEFAULTTYPE_API BoundingBox1D
{

public:
    typedef std::pair< SReal, SReal > bbox_t;

    BoundingBox1D();
    /// Define using the endpoints of the main diagonal
    BoundingBox1D(SReal minBBox, SReal maxBBox);
    BoundingBox1D(const bbox_t& bbox);
    /// Define using xmin, xmax in this order
    BoundingBox1D(const Vec2f& bbox);
    /// Define using xmin, xmax in this order
    BoundingBox1D(const Vec2d& bbox);

    static BoundingBox1D neutral_bbox();

    operator bbox_t() const;

    void invalidate();
    bool isValid() const;
    bool isFlat()  const;
    bool isNegligeable() const; // !valid || flat
    bool isNull()  const;

    const SReal&  minBBox() const;
    const SReal&  maxBBox() const;
    SReal& minBBox();
    SReal& maxBBox();

    bool contains( SReal point) const;
    bool contains( const BoundingBox1D& other) const;

    bool intersect( const BoundingBox1D& other) const;
    void intersection( const BoundingBox1D& other);

    void include( SReal point);
    void include( const BoundingBox1D& other);

    void inflate( SReal amount);

    BoundingBox1D getIntersection( const BoundingBox1D& other ) const;
    BoundingBox1D getInclude( SReal point ) const;
    BoundingBox1D getInclude( const BoundingBox1D& other ) const;
    BoundingBox1D getInflate( SReal amount ) const;

    friend std::ostream& operator << ( std::ostream& out, const BoundingBox1D& bbox)
    {
        out << bbox.minBBox() << " " <<  bbox.maxBBox();
        return out;
    }

    friend std::istream& operator >> ( std::istream& in, BoundingBox1D& bbox)
    {
        in >> bbox.minBBox() >> bbox.maxBBox();
        return in;
    }


protected:
    bbox_t bbox;
};

template <>
struct DataTypeName<BoundingBox>
{
  static const char *name() { return "BoundingBox"; }
};

struct BoundingBoxTypeInfo
{
    typedef BoundingBox DataType;
    typedef Vector3 BaseType;
    typedef DataTypeInfo<Vec3d> BaseTypeInfo;
    typedef typename BaseTypeInfo::ValueType ValueType;

    enum
    {
        ValidInfo = BaseTypeInfo::ValidInfo
    };  ///< 1 if this type has valid infos
    enum
    {
        FixedSize = 1
    };  ///< 1 if this type has a fixed size  -> always 1 single pair of vec3
    enum
    {
        ZeroConstructor = 0
    };  ///< 1 if the constructor is equivalent to setting memory to 0
    /// -> I don't think so, bbox are initialized with +inf / -inf usually
    enum
    {
        SimpleCopy = 1
    };  ///< 1 if copying the data can be done with a memcpy
    enum
    {
        SimpleLayout = 1
    };  ///< 1 if the layout in memory is simply N values of the same base type
    enum
    {
        Integer = 0
    };  ///< 1 if this type uses integer values
    enum
    {
        Scalar = 1
    };  ///< 1 if this type uses scalar values
    enum
    {
        Text = 0
    };  ///< 1 if this type uses text values
    enum
    {
        CopyOnWrite = 1
    };  ///< 1 if this type uses copy-on-write -> it seems to be THE important
    ///< option not to perform too many copies
    enum
    {
        Container = 1
    };  ///< 1 if this type is a container

    enum
    {
        Size = 2
    };  ///< largest known fixed size for this type, as returned by size()

    static size_t size() { return 3; } // supposed to be the total number of elements. Ends up being the number of elements in the 2nd dimension
    static size_t byteSize() { return sizeof (ValueType); }  // Size of the smalest single element in the container: BoundingBox uses Vec3d internally, so double

    static size_t size(const DataType & /*data*/) { return 2 * BaseTypeInfo::size(); } // supposed to be the nb of elements in the 1st dimension. Ends up being the total number of elems.

    static bool setSize(DataType & /*data*/, size_t /*size*/) { return false; } // FixedArray -> ignore

    template <typename T>
    static void getValue(const DataType & data, size_t index,
                         T & value) /// since TypeInfos abstract all containers as 1D arrays, T here is of ValueType
    {
        value = static_cast<T>(((ValueType*)&data)[index]);
    }

    template <typename T>
    static void setValue(DataType & data, size_t index,
                         const T & value)
    {
        ((ValueType*)&data)[index] = static_cast<ValueType>(value);
    }

    static double getScalarValue (const void* data, size_t index)
    {
        return ((ValueType*)&data)[index];
    }

    static void setScalarValue (const void* data, size_t index, double value)
    {
        ((ValueType*)&data)[index] = value;
    }


    static void getValueString(const DataType &data, size_t index,
                               std::string &value)
    {
        if (index != 0) return;
        std::ostringstream o;
        o << data;
        value = o.str();
    }

    static void setValueString(DataType &data, size_t index,
                               const std::string &value)
    {
        if (index != 0) return;
        std::istringstream i(value);
        i >> data;
    }

    static const void *getValuePtr(const DataType & bbox) { return (const void*)(&bbox); }

    static void *getValuePtr(DataType &bbox) { return (void*)(&bbox); }
};

template <>
struct DataTypeInfo<BoundingBox> : public BoundingBoxTypeInfo
{
    static std::string name() { return "BoundingBox"; }
};


}
}

#endif // SOFA_DEFAULTTYPE_BOUNDINGBOX_H
