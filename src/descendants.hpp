//
// 
// Code to assign descendant information to event.
// 
//

#ifndef DESCENDANTS_H
#define DESCENDANTS_H


#include "scrm/src/model.h"



typedef uint64_t Descendants_t;

static const Descendants_t NO_DESCENDANTS = 0;



inline Descendants_t get_descendants( const Node* node ) {

    Descendants_t result = 0;
    while (node->label() == 0) {
        result |= get_descendants( node->getLocalChild1() );
        node = node->getLocalChild2();  // could be NULL
        if (!node) return result;
    }
    assert (node->label() <= 64);
    return result | (((Descendants_t)1) << (node->label() - 1));
}


// updates 'sample' to next 1-based descendant index, or returns 'false'
// initial call should be with sample == 0
inline bool get_next_descendant( Descendants_t& descendants, int& sample ) {

    if (descendants == 0) return false;
    while ( (descendants & 1) == 0 ) {
        ++sample;
        descendants >>= 1;
    }
    ++sample;
    descendants >>= 1;
    return true;
}


#endif
