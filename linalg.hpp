#include <math.h>
#include <initializer_list>

namespace linalg {

    template <typename T>
    class vec3;

    template <typename T>
    class ray3;

    template <typename T>
    class plane3;

    template <typename T, unsigned int N>
    class poly;    

    //  Simple 3D vector class
    template <typename T>
    class vec3 {
    public:
        T x;
        T y;
        T z;
        //T norm;

        vec3(T x, T y, T z): x(x), y(y), z(z) {}

        vec3(const T* coords): x(*coords), y(*(coords + 1)), z(*(coords + 2)) {}
        
        vec3(const vec3& other) {
            *this = other;
        }

        vec3(): x(0), y(0), z(0) {}

        inline T norm() {
            return sqrt(x * x + y * y + z * z);
        }

        inline T cos(const vec3& other) {
            return *this * other / norm() / other.norm();
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
        }

        inline void normalize() {
            *this /= norm();
        }

        inline bool operator == (const vec3& other) {
            return this->x == other.x && this->y == other.y && this->z == other.z;
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
            u.normalize();
        }

        ray3(T p0x, T p0y, T p0z, T ux, T uy, T uz): u(ux, uy, uz), p0(p0x, p0y, p0z) {
            u.normalize();
        }

        ray3(vec3<T> &p0, vec3<T> &u): u(u), p0(p0) {
            u.normalize();
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

        template <unsigned int N>
        T intersect_poly(poly<T, N>& p) {
            return N;
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
            n.normalize();
        }

        plane3(T v0x, T v0y, T v0z, T norm_x, T norm_y, T norm_z): v0(v0x, v0y, v0z), n(norm_x, norm_y, norm_z) {
            n.normalize();
        }

        plane3(vec3<T> &v0, vec3<T> &norm): n(norm), v0(v0) {
            n.normalize();
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

    template <typename T, unsigned int N>
    class poly {
    private:
        void adder(T* arr, int idx, size_t sz, T v) {
            arr[idx] = v;
        }

        template <typename...Args>
        void adder(T* arr, int idx, size_t sz, T fst, Args... args) {
            arr[idx++] = fst;
            if (idx < sz) {
                adder(arr, idx, sz, args...);
            }
        }

    public:        
        vec3<T> points[N];
        vec3<T> norm;

        template <typename...Args>
        poly(Args... args) {
            auto sz = N * sizeof(vec3<T>) / sizeof(T);
            T* loc_points = (T*)points;
            adder(loc_points, 0, sz, args...);

            for (int i = 0; i < N; ++i) {
                printf("%f, %f, %f\n", points[i].x, points[i].y, points[i].z);
            }
        }

        poly(T src_points[]) {
            memcpy(points, src_points, sizeof(T) * N);
            norm = (points[1] - points[0]).cross(points[N] - points[0]);
        }

        T get_area() {
            T result = 0;
            vec3<T> &zero_p = points[0];
            for (int i = 1; i < N - 1; ++i) {
                vec3<T> a = points[i] - zero_p;
                vec3<T> b = points[i + 1] - zero_p;
                result += a.cross(b) * 0.5;
            }
            return result;
        }

        bool check_plane() {
            vec3<T> &zero_p = points[0];
            vec3<T> cross = (points[1] - zero_p).cross(points[2] - zero_p);
            for (int i = 2; i < N - 1; ++i) {
                vec3<T> a = points[i] - zero_p;
                vec3<T> b = points[i + 1] - zero_p;
                bool eq = a.cross(b) == cross;
                if (!eq) {
                    return false;
                }
            }
            return true;
        }
    };

}