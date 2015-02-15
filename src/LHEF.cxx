#include "LHEF.h"
// Les Houches Event File class.
using namespace LHEF;

// Destructor implemented here.

void EventGroup::clear() {
  while ( size() > 0 ) {
    delete back();
    pop_back();
  }
};

EventGroup::~EventGroup() {
  clear();
};

EventGroup::EventGroup(const EventGroup &) {
  for ( int i = 0, N = size(); i < N; ++i ) at(i) = new HEPEUP(*at(i));
};

EventGroup & EventGroup::operator=(const EventGroup & x) {
  if ( &x == this ) return *this;
  clear();
  nreal = x.nreal;
  ncounter = x.ncounter;
  for ( int i = 0, N = x.size(); i < N; ++i ) push_back(new HEPEUP(*x.at(i)));
  return *this;
};
