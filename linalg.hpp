#include <math.h>

namespace linalg {

    template <typename T>
    class vec3;

    template <typename T>
    class ray3;

    template <typename T>
    class plane3;

    //  Simple 3D vector class
    template <typename T>
    class vec3 {
    public:
        T x;
        T y;
        T z;
        T norm;

        vec3(T x, T y, T z): x(x), y(y), z(z), norm(sqrt(x * x + y * y + z * z)) {}

        vec3(const T* coords): x(*coords), y(*(coords + 1)), z(*(coords + 2)) {
            norm = sqrt(x * x + y * y + z * z);
        }
        
        vec3(const vec3& other) {
            *this = other;
        }

        vec3(): x(0), y(0), z(0), norm(0) {}

        inline T cos(const vec3& other) {
            return *this * other / norm / other.norm;
        }

        inline vec3 cross(const vec3& other) {
            return vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
        }

        inline vec3 operator + (const vec3& other) {
            return vec3(x + other.x, y + other.y, z + other.z);
        }

        inline vec3 operator - (const vec3& other) {
            return vec3(x - other.x, y - other.y, z - other.z);
        }

        inline T operator * (const vec3& other) {
            return x * other.x + y * other.y + z * other.z;
        }

        inline vec3 operator * (T a) {
            return vec3(a * x, a * y, a * z);
        }

        inline vec3 operator / (T a) {
            return vec3(x / a, y / a, z / a);
        }

        inline void operator /= (T a) {
            x /= a;
            y /= a;
            z /= a;
            norm /= a;
        }
    };

    // 3D ray class
    template <typename T>
    class ray3 {
    public:
        // Direction
        vec3<T> u;
        // Start-point
        vec3<T> p0;

        ray3(const T* coords): p0(coords), u(coords + 3) {
            u /= u.norm;
        }

        ray3(T p0x, T p0y, T p0z, T ux, T uy, T uz): u(ux, uy, uz), p0(p0x, p0y, p0z) {
            u /= u.norm;
        }

        ray3(vec3<T> &p0, vec3<T> &u): u(u), p0(p0) {
            u /= u.norm;
        }

        T intersect_plane(plane3<T>& plane, vec3<T>& intersection) {
            T denom = plane.n * u;
            if (denom >= 0) {
                return -1;
            }
            T s = plane.n * (plane.v0 - p0) / denom;
            intersection = p0 + u * s;
            return s;
        }
    };

    // 3D plane (hyperplane) class
    template <typename T>
    class plane3 {
    public:    
        // normal
        vec3<T> n;
        // Norm-point
        vec3<T> v0;

        plane3(const T* coords): v0(coords), n(coords + 3) {
            n /= n.norm;
        }

        plane3(T v0x, T v0y, T v0z, T norm_x, T norm_y, T norm_z): v0(v0x, v0y, v0z), n(norm_x, norm_y, norm_z) {
            n /= n.norm;
        }

        plane3(vec3<T> &v0, vec3<T> &norm): n(norm), v0(v0) {
            n /= n.norm;
        }

        T intersect_ray(ray3<T>& ray, vec3<T>& intersection) {
            T denom = n * ray.u;
            if (denom >= 0) {
                return -1;
            }
            T s = n * (v0 - ray.p0) / denom;
            intersection = ray.p0 + ray.u * s;
            return s;
        }
    };

}