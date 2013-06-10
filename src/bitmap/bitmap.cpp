#include "bitmap.h"

void BitmapIndex::loadIndex(const char* fname) {
    
}
void BitmapIndex::saveIndex(const char* fname) {
    
}

template <class T>
BitSlicedIndex<T>::BitSlicedIndex(unsigned int nwords, const char* fname) {
    
}

template <class T>
void BitSlicedIndex<T>::append(T value) {
    
}
template <class T>
void BitSlicedIndex<T>::append(T* value) {
    
}
template <class T>
T BitSlicedIndex<T>::decode(size_t idx) {
    
}
template <class T>
T* BitSlicedIndex<T>::decode(size_t idx) {
    
}

template <class T>
RangeEncodedIndex<T>::RangeEncodedIndex(vector<T> vec, const char* fname) {
    
}

template <class T>
void RangeEncodedIndex<T>::append(T value) {
    
}

template <class T>
T RangeEncodedIndex<T>::decode(size_t idx) {
    
}
