/* 
    This file is a part of eXlibris C++ Library
    under the GNU General Public License:
    See the LICENSE.md files for terms and 
    conditions.
*/


#ifndef _SUBGROUPPARTMAN_H
#define _SUBGROUPPARTMAN_H

// std
#include <fstream>
#include <cassert>
#include <cmath>
// xtools
#include "xMPIEnv.h"
#include "xPartitionManagerTools.h"
//xfem::distmesh
#include "xLoadBalanceTools.h"
//xmeshtool
#include "xSplitAOMD.h"
// xexport
#include "xExportGmsh.h"
#include "xExportAlgorithm.h"
// xfem
#include "xData.h"
#include "xAlgorithm.h"
#include "xExportGmsh.h"


using namespace xtool;
using namespace xfem;
using namespace std;
using namespace AOMD;
using namespace xexport;

//#define RANDOM

template < typename T >
using datamanager_t = xinterface::aomd::xAttachedDataManagerAOMD < T >;

// Level criteria to get uniform refinement
// => transfert+criteria
class transLevel : public xinterface::xtemplaterefinemesh::xSplitTransferBaseAOMD
{
    public:
        transLevel(unsigned int tag_level_) : tag_level(tag_level_) {}
        void collect(AOMD::mEntity *e)
        {
            if (e->getLevel() == 3)
                cur_level = e->getAttachedInt(tag_level);
        }
        void transfer(AOMD::mEntity *se)
        {
            if (se->getLevel() == 3)
                se->attachInt(tag_level,cur_level+1);
        }
    private:
        int cur_level;
        unsigned int tag_level;
};


class criteriaLevel
{
    public:
        criteriaLevel(unsigned int tag_level_,int max_level_) :
            tag_level(tag_level_)
            ,max_level(max_level_)
        {}
        void update ( void)
        {}
        bool operator () (AOMD::mEntity &e)
        {
            if ( e.getAttachedInt(tag_level) < max_level ) return true;
            else return false;
        }
    private:
        unsigned int tag_level;
        int max_level;
};
class SGIDInfoManager
{
    public:
        SGIDInfoManager(const datamanager_t < int > & sub_group_,datamanager_t < bool > & sub_group_remote_) :
            sub_group(sub_group_)
            ,sub_group_remote(sub_group_remote_)
        {}
        // mandatory traits
        typedef  keyManagerSendAndReceive::information_key_t information_key_t;
        typedef unsigned short information_t;
        typedef xtool::homogeneous_data_style_trait data_style_trait;
        typedef xtool::send_and_recv_keys_communication_trait communication_trait;

        // mandatory methods
        void setInfo(information_key_t key, const information_t &info, int receivedfrom)
        {
            AOMD::mEntity *f=const_cast<AOMD::mEntity *>(key);
            assert(key->size(3) == 1);
            AOMD::mEntity * tet = key->get(3,0);
            const int *color = sub_group.getData(*tet);
            assert(color);
            if ( *color == info)
            {
                sub_group_remote.setData(*f) = true;
            }
            else
            {
                sub_group_remote.deleteData(*f);
            }
        }
        information_t getInfo(information_key_t key, int sendto)
        {
            assert(key->size(3) == 1);
            AOMD::mEntity * tet = key->get(3,0);
            const int *color = sub_group.getData(*tet);
            assert(color);
            return *color;
        }
    private:
        const datamanager_t < int > & sub_group;
        datamanager_t < bool > & sub_group_remote;

};


#endif



