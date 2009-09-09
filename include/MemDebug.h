///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
#ifndef MemDebugH
#define MemDebugH
//---------------------------------------------------------------------------


/* Put the following in the main program file:

#ifdef MEMDEBUG
AllocList* allocList=0;
#endif

and somewhere in main:

atexit(&DumpUnfreed);

*/

#ifdef MEMDEBUG

// Here is some code for tracking memory leaks I got from
// http://www.flipcode.com/tutorials/tut_memleak.shtml
// Code is originally (and still mostly) by Dion Picco (c) 23 May 2000

// Remember to use atexit(&DumpUnfreed) to print results

#include <string>
#include <iostream>

struct AllocInfo {
  AllocInfo(long s, const char* f, int l) :
      size(s), file(f), line(l) {}

  long size;
  std::string file;
  int line;
};

#include <map>

typedef std::map<void*,AllocInfo*> AllocList;

extern AllocList *allocList;

inline void AddTrack(void* addr, long size, const char* file, int line)
{
  if(!allocList) {
    allocList = new AllocList;
  }

  if (allocList->count(addr)) {
    std::cout<<addr<<" already in allocList, trying to add:\n";
    std::cout<<"size = "<<size<<std::endl;
    std::cout<<"file = "<<file<<", line num = "<<line<<std::endl;
    abort();
  }
  (*allocList)[addr] = new AllocInfo(size,file,line);
}

//inline void RemoveTrack(void* addr, const char* file, int line)
inline void RemoveTrack(void* addr)
{
  if(!allocList || !allocList->count(addr)) {
    //std::cout<<"deleting "<<addr<<", but not in AllocList\n";
    //std::cout<<"file = "<<file<<", line num = "<<line<<std::endl;
    return;
  }
  allocList->erase(addr);
}

inline void DumpUnfreed()
{
  long totalSize = 0;
  if(!allocList) return;
  for(AllocList::iterator i=allocList->begin();i!=allocList->end();++i) {
    std::cout << i->second->file <<"  LINE "<<i->second->line<<",\t\tADDRESS ";
    std::cout << i->first <<"  "<<i->second->size<<" bytes unfreed\n";
    totalSize += i->second->size;
  }
  std::cout << "-----------------------------------------------------------\n";
  std::cout << "Total Unfreed: "<<totalSize<<" bytes\n";
}

inline void* operator new(size_t size, const char* file,int line)
{ 
  void* ptr = (void*)malloc(size);
  AddTrack(ptr, size, file, line);
  return(ptr);
}

//inline void operator delete(void* p, const char* file, int line)
void operator delete(void* p) throw()
{
  //RemoveTrack(p,file,line);
  RemoveTrack(p);
  free(p);
}

inline void* operator new[](size_t size,const char* file,int line)
{ 
  void* ptr = (void*)malloc(size);
  AddTrack(ptr, size, file, line);
  return(ptr);
}

//inline void operator delete[](void* p, const char* file, int line)
void operator delete[](void* p) throw()
{
  //RemoveTrack(p,file,line);
  RemoveTrack(p);
  free(p);
}


#define MEMDEBUG_NEW new(__FILE__, __LINE__)
#define new MEMDEBUG_NEW

//#define MEMDEBUG_DELETE delete(__FILE__, __LINE__)
//#define delete MEMDEBUG_DELETE

#endif // ifdef MEMDEBUG

#endif // MemDebug_H
