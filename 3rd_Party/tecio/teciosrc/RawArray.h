 #pragma once
#include <algorithm>
#include "ClassMacros.h"
#include "CodeContract.h"
#include "MiscMacros.h"
namespace tecplot { template <typename T> class ___3269 { public: void ___2319( char const* container, size_t numElements) {
 #ifdef LARGE_ARRAY_MEMORY_LOGGING
size_t const MEMTRACK_CUTOFF = size_t(1000)*size_t(1000); if ( numElements * sizeof(T) >= MEMTRACK_CUTOFF ) { FILE *file = fopen("memtrack.txt", "at"); if ( file ) { fprintf(file, "%s\t%" "I64u" "\t%" "I64u" "\t%s\n", container, numElements, sizeof(T), typeid(T).___2685()); fclose(file); } else throw std::bad_alloc(); }
 #else
___4278(container); ___4278(numElements);
 #endif
} typedef T value_type;
 #if (defined _MSC_VER && __cplusplus >= 199711L) || __cplusplus >= 201103L
___3269(___3269 const&) = delete; ___3269& operator=(___3269 const&) = delete;
 #if defined _MSC_VER && _MSC_VER <= 1800 
___3269(___3269&& ___2888) : m_valuesRef(std::move(___2888.m_valuesRef)), ___2671(std::move(___2888.___2671)), m_capacity(std::move(___2888.m_capacity)), m_size(std::move(___2888.m_size)) {} ___3269& operator=(___3269&& ___3392) { if (this != &___3392) { m_valuesRef = std::move(___3392.m_valuesRef); ___2671 = std::move(___3392.___2671); m_capacity = std::move(___3392.m_capacity); m_size = std::move(___3392.m_size); } return *this; }
 #else
___3269(___3269&&) = default; ___3269& operator=(___3269&&) = default;
 #endif
 #else
private: ___3269(___3269 const&); ___3269& operator=(___3269 const&); public:
 #endif
inline ___3269() : m_valuesRef(NULL) , ___2671(NULL) , m_capacity(0) , m_size(0) { } inline ___3269( T* const& values, size_t    capacity, size_t    size = 0) : m_valuesRef(NULL) , ___2671(const_cast<T*>(values)) , m_capacity(capacity) , m_size(size) { REQUIRE(VALID_REF(this->___2671)); REQUIRE(m_capacity != 0); REQUIRE(this->m_capacity >= this->m_size); } inline ___3269( T*&    values, size_t capacity = 0, size_t size = 0) : m_valuesRef(&values) , ___2671(values) , m_capacity(capacity) , m_size(size) { REQUIRE(VALID_REF(this->___2671) || this->___2671 == NULL); REQUIRE(EQUIVALENCE(this->___2671 == NULL, m_capacity == 0)); REQUIRE(this->m_capacity >= this->m_size); } inline void nullify() { if (m_valuesRef != NULL) *m_valuesRef = NULL; m_valuesRef = NULL; ___2671 = NULL; m_capacity = 0; m_size = 0; } inline bool empty() const { return m_size == 0; } inline size_t size() const { return m_size; } inline void reserve(size_t size) { if (size > m_capacity) enlargeCapacity(size); } inline void ___3503(size_t size) { REQUIRE(size <= m_capacity); m_size = size; ENSURE(m_size <= m_capacity); } inline void clear() { m_size = 0; } inline size_t capacity() const { return m_capacity; } inline ___3269& copy(___3269 const& ___2888) { reserve(___2888.m_size); ___3503(___2888.m_size); for (size_t ___2865 = 0; ___2865 < m_size; ++___2865) ___2671[___2865] = ___2888.___2671[___2865]; return *this; } inline ___3269& append(___3269 const& ___2888) { size_t const origSize = m_size; reserve(m_size + ___2888.m_size); ___3503(m_size + ___2888.m_size); for (size_t ___2865 = origSize; ___2865 < m_size; ++___2865) ___2671[___2865] = ___2888.___2671[___2865-origSize]; return *this; } inline ___3269& append(T const& ___4314) { size_t const origSize = m_size; reserve(m_size + 1); ___3503(m_size + 1); ___2671[origSize] = ___4314; return *this; } inline void copyTo(T* target) { for (size_t ___2865 = 0; ___2865 < m_size; ++___2865) target[___2865] = ___2671[___2865]; } inline T& operator[](size_t ___2865) { REQUIRE(___2865 < this->m_size || (___2865 == 0 && this->m_capacity != 0)); return ___2671[___2865]; } inline T const& operator[](size_t ___2865) const { REQUIRE(___2865 < this->m_size || (___2865 == 0 && this->m_capacity != 0)); return ___2671[___2865]; } inline bool ___2035() const { return (m_valuesRef == NULL || *m_valuesRef == NULL); } typedef T const* const_iterator; inline const_iterator begin() const { return ___2671; } inline const_iterator end() const { return ___2671+m_size; } typedef T* iterator; inline iterator begin() { return ___2671; } inline iterator end() { return ___2671+m_size; } private: inline void enlargeCapacity(size_t capacity) { REQUIRE(capacity > m_capacity); REQUIRE(m_valuesRef != NULL); ___2319("RawArray", capacity); T* values = new T[capacity]; if (m_size != 0) std::copy(___2671, ___2671+m_size, values); delete [] ___2671; *m_valuesRef = values; ___2671     = values; m_capacity   = capacity; } T**    m_valuesRef; T*     ___2671; size_t m_capacity; size_t m_size; }; }
