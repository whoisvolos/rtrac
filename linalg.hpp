#include <math.h>
#include <tuple>

namespace linalg {
    float EPS = 1e-7;

    template <typename T>
    class vec3;

    template <typename T>
    class ray3;

    template <typename T>
    class plane3;

    template <typename T, size_t N>
    class poly;

    template <typename T>
    class vec2 {
    public:
        T u;
        T v;

        vec2(): u(0), v(0) {}

        vec2(T u, T v): u(u), v(v) {}

        inline vec2 operator + (const vec2& other) {
            return vec2(u + other.u, v + other.v);
        }

        inline vec2 operator - (const vec2& other) {
            return vec2(u - other.u, v - other.v);
        }

        inline T norm() {
            return sqrt(u * u + v * v);
        }

        inline T perp(const vec2& other) {
            return u * other.v - other.u * v;
        }

        inline static vec2<T> from3(const vec3<T>& vec, int drop_idx) {
            switch (drop_idx) {
                case 0:
                    return vec2(vec.y, vec.z);
                case 1:
                    return vec2(vec.x, vec.z);
                default:
                case 2:
                    return vec2(vec.x, vec.y);
            }
        }
    };

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
            return fabs(this->x - other.x) < EPS &&
                   fabs(this->y - other.y) < EPS &&
                   fabs(this->z - other.z) < EPS;
        }

        inline std::tuple<int, T> max() {
            int idx = max_magn_idx();
            return std::make_tuple(idx, *((T *)this + idx));
        }

        inline int max_magn_idx() {
            T xa = abs(x), ya = abs(y), za = abs(z);
            T xy = xa > ya ? xa : ya;
            int idx = xa > ya ? 0 : 1;
            return xy > za ? idx : 2;
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
            if (fabs(denom) < EPS) {
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
            if (fabs(denom) < EPS) {
                return -1;
            }
            T s = n * (v0 - ray.p0) / denom;
            intersection = ray.p0 + ray.u * s;
            return s;
        }
    };

    template <typename T, size_t N>
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

        poly(T* src_points) {
            memcpy(points, src_points, sizeof(vec3<T>) * N);
        }

        template <typename...Args>
        poly(Args... args) {
            auto sz = N * sizeof(vec3<T>) / sizeof(T);
            T* loc_points = (T*)points;
            adder(loc_points, 0, sz, args...);
        }

        vec3<T> norm() {
            return (points[1] - points[0]).cross(points[N - 1] - points[0]);
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

        T intersect_ray(ray3<T>& ray, vec3<T>& intersection) {
            // Find plane intersection
            vec3<T> norm = this->norm();
            plane3<T> plane(points[0], norm);
            T dist = ray.intersect_plane(plane, intersection);
            if (fabs(dist) < EPS) {
                return 0;
            }

            // Drop maximum magnitude coord
            auto max_idx = norm.max_magn_idx();
            vec2<T> projection[N];
            vec2<T> proj_intr = vec2<T>::from3(intersection, max_idx);
            //printf("After %i drop: %f, %f\n", max_idx, proj_intr.u, proj_intr.v);
            for (int i = 0; i < N; ++i) {
                projection[i] = vec2<T>::from3(*(points + i), max_idx) - proj_intr;
                //printf("%f, %f\n", projection[i].u, projection[i].v);
            }

            int intrs = 0;
            vec2<T> u(1, 0);
            for (int i = 0; i < N; ++i) {
                vec2<T> w = projection[i];
                vec2<T> v = projection[(i + 1) % N] - w;
                T denom = u.perp(v);
                if (fabs(denom) >= EPS) {
                    T numer = u.perp(w);
                    printf("%f / %f\n", numer, denom);
                    if (numer / denom <= v.norm()) {
                        ++intrs;
                    }
                }
            }

            if (intrs % 2 == 0) {
                printf("Intersection!\n");
            }

            return dist;
        }
    };

}