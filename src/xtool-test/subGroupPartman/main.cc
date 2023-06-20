/*
    This file is a part of eXlibris C++ Library
    under the GNU General Public License:
    See the LICENSE.md files for terms and
    conditions.
 */

#include "main.h"

#include "xMesh.h"
int main(int argc, char *argv[])
{
   xMPIEnv::init(argc, argv);

   int proc_id, nb_proc, nb_sg;
   char names[256];
   MPI_Comm world = MPI_COMM_WORLD;
   MPI_Comm_rank(world, &proc_id);
   MPI_Comm_size(world, &nb_proc);
#ifndef RANDOM
   if (nb_proc != 4)
   {
      cout << "Test designed to run on 4 process in non RAMDOM mode !" << endl;
      throw -476;
   }
#endif
   nb_sg = 2 * nb_proc;
   xEvalConstant<double> eval_pid(proc_id);

   // Redirection of all output to individual files
   string pid = std::to_string(proc_id);
   string fo = "proc_" + pid + "_output.txt";
   FILE *ok = freopen(fo.c_str(), "w", stdout);
   if (!ok)
   {
      std::cout << "Can't reopen stdout on file " << fo << __FILE__ << __LINE__ << " compiled " << __DATE__ << " at " << __TIME__
                << std::endl;
      throw;
   }

   cout << "Starting reading the master data file" << endl;
   xfem::xData data;
   data.ReadInfo("data/main.dat");
   cout << "Done with reading the master data file" << endl;

   cout << "Starting reading the mesh file" << endl;
   datamanager_t<int> target_proc;
   data.ReadMesh(target_proc, world);
   if (!proc_id && target_proc.beginKey() == target_proc.endKey())
   {
      cout << "No available partition !" << endl;
      throw -456;
   }
   cout << "End reading the mesh file" << endl;

   xMesh &m = *data.getMesh();
   auto &m_pm = m.getPartitionManager();

   // clean before migration
   data.entities_id.clear();

   // dispatch original mesh
   cout << "Starting coarse migration" << endl;
   xmeshtool::migrateEntities(m.getMesh(), m_pm, target_proc);
   target_proc.clear();
   cout << "End coarse migration" << endl;

#ifndef RANDOM
   xgeom::xBoundingBox bb = m.compute_bounding_box();
#endif

   cout << "Starting fine mesh creation" << endl;
   // set Level Criteria
   unsigned int tag_level = AOMD::AOMD_Util::Instance()->newMeshDataId("mytaglevel");
   criteriaLevel level_criteria(tag_level, 4);
   transLevel trans(tag_level);

   xinterface::xtemplaterefinemesh::xMeshSplitUserCriteriaAOMD<criteriaLevel> splitor(level_criteria, 6);
   splitor.registerTransfer(&trans);

   // split mesh
   splitor.split(m.getMesh(), 3, m_pm);
   cout << "End fine mesh creation" << endl;

#ifdef RANDOM
   cout << "Starting fine partition creation" << endl;
   // set  fine partition with parmetis
   datamanager_t<int> fine_target_proc = xmeshtool::setParmetisPartition(m, m_pm, nb_proc);
   cout << "End fine partition creation" << endl;

   // dispatch fine mesh
   cout << "Starting fine migration" << endl;
   xmeshtool::migrateEntities(m, m_pm, fine_target_proc);
   fine_target_proc.clear();
   cout << "End fine migration" << endl;
#endif

   // export working mesh
   string expfine = "fine_partion_" + std::to_string(proc_id) + ".msh";
   AOMD::AOMD_Util::Instance()->exportGmshFile(expfine.c_str(), &m.getMesh());

   cout << "Starting sub group creation" << endl;
#ifdef RANDOM
   // use parmetis to set groups
   datamanager_t<int> sub_group = xmeshtool::setParmetisPartition(m, m_pm, nb_sg);
#else
   datamanager_t<int> sub_group;
   /*
    * oblique plane
      double a = bb.max(0)-bb.min(0);
      double b = bb.max(1)-bb.min(1);
      double d = -a *bb.min(0)-b*bb.min(1);
      double r = fabs(a*bb.max(0)+b*bb.max(1)+d);
      if (r == 0.0)
      {
       cout<<"bad geom"<<endl;
       throw -23.;
      }
      r = 1./r;
      std::for_each(m.begin(3),m.end(3),[&a,&b,&d,&r,&nb_sg,&m,&sub_group](mEntity *e)
                 {
                     xtensor::xPoint cg = m.COG(e);
                     double dp=fabs(a*cg(0)+b*cg(1)+d)*r;
                     sub_group.setData(*e)=(static_cast<int>(floor(dp*nb_sg)))%nb_sg;
                 });
    */
   // Sphere pattern for groups. This setting is giving best variety:
   //    - edges are present in group partition managers
   //    - element are present in group only in one process at general process interface => filtering
   //    - all groups got a partition manager
   double r = bb.max(0) * bb.max(0) + bb.max(1) * bb.max(1);
   if (r == 0.0)
   {
      cout << "bad geom" << endl;
      throw -23.;
   }
   r = 1. / r;
   std::for_each(m.begin(3), m.end(3), [&r, &nb_sg, &m, &sub_group](mEntity *e) {
      Trellis_Util::mPoint cg = m.getMesh().COG(e);
      double dp = (cg(0) * cg(0) + cg(1) * cg(1) + cg(2) * cg(2)) * r;
      sub_group.setData(*e) = (static_cast<int>(floor(dp * nb_sg))) % nb_sg;
   });
#endif

   // export them
   xEvalFromLambda<double> evalSG([&sub_group](const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ) {
      AOMD::mEntity *e_appro = geo_appro->getEntity();
      int *color = sub_group.getData(*e_appro);
      assert(color);
      return *color;
   });
   // xexport::xExportGmshMPIIO mexport(world);
   // xexport::xExportGmshAscii mexport(world);
   xexport::xExportGmshAsciiSort mexport(world);
   sprintf(names, "sub_group");
   xIntegrationRuleBasic integrule_basic(0);
   Export(evalSG, mexport, names, integrule_basic, m.begin(3), m.end(3));
   cout << "End sub group creation" << endl;

   cout << "Starting face boundary connected to element of same group creation" << endl;
   // preparation: exchange sg_id on face proc boundary
   // create key container
   keyManagerSendAndReceive keymanager(m_pm);
   xKeyContainerSendAndRecv<keyManagerSendAndReceive::information_key_t> keycontainer(world);

   // Set the key: only interested by faces
   // use sub_group_remote as intermediate
   datamanager_t<bool> sub_group_remote;
   vector<keyManagerSendAndReceive::information_key_t> faces;
   for (auto e : make_range(m_pm.beginObject(), m_pm.endObject()))
   {
      if (e->getLevel() == 2) sub_group_remote.setData(*e) = true;
   }
   keycontainer.accumulateKeysAllGather(sub_group_remote.beginKey(), sub_group_remote.endKey(), keymanager);

   // Exchange for all boundary faces their SG ID counterpart and check with local: store in sub_group_remote
   // only face having same SG ID counterpart.
   SGIDInfoManager infomanager(sub_group, sub_group_remote);
   exchangeInformation(keycontainer, infomanager);
   cout << "Ending face boundary connected to element of same group creation" << endl;

   cout << "Starting test function createPartitionManagerForSubGroup" << endl;
   // for first nb_proc sub group we are going to consider the original communicator
   int i;
   for (i = 0; i < nb_proc; ++i)
   {
      partmanAOMD_t part_man_SG(world);
      xFilteredRegion<datamanager_t<bool>::c_iterKey_t, xEntityFilter> fr(sub_group_remote.beginKey(), sub_group_remote.endKey(),
                                                                          [&sub_group, &i](AOMD::mEntity *f) -> bool {
                                                                             assert(f->getLevel() == 2);
                                                                             assert(f->size(3) == 1);
                                                                             AOMD::mEntity *tet = f->get(3, 0);
                                                                             int *color = sub_group.getData(*tet);
                                                                             assert(color);
                                                                             if (*color == i) return true;
                                                                             return false;
                                                                          });
      createPartitionManagerForSubGroup<datamanager_t, AOMD::mEntity>(make_range(fr.begin(), fr.end()), m_pm, part_man_SG);
      sprintf(names, "sub_group_part_man_%d", i);
      Export(eval_pid, mexport, names, integrule_basic, part_man_SG.beginObject(), part_man_SG.endObject());
   }
   // for last set of group we are going to consider a restricted communicator: the odd/even communicator
   MPI_Comm odd_even_world = MPI_COMM_WORLD;
   MPI_Comm_split(world, proc_id % 2, 0, &odd_even_world);
   xexport::xExportGmshAsciiSort oemexport(odd_even_world);
   for (; i < nb_sg; ++i)
   {
      partmanAOMD_t part_man_SG(odd_even_world);
      xFilteredRegion<datamanager_t<bool>::c_iterKey_t, xEntityFilter> fr(sub_group_remote.beginKey(), sub_group_remote.endKey(),
                                                                          [&sub_group, &i](AOMD::mEntity *f) -> bool {
                                                                             assert(f->getLevel() == 2);
                                                                             assert(f->size(3) == 1);
                                                                             AOMD::mEntity *tet = f->get(3, 0);
                                                                             int *color = sub_group.getData(*tet);
                                                                             assert(color);
                                                                             if (*color == i) return true;
                                                                             return false;
                                                                          });
      createPartitionManagerForSubGroup<datamanager_t, AOMD::mEntity>(make_range(fr.begin(), fr.end()), m_pm, part_man_SG);
      sprintf(names, "sub_group_part_man_%d_%d", proc_id, i);
      Export(eval_pid, oemexport, names, integrule_basic, part_man_SG.beginObject(), part_man_SG.endObject());
   }
   cout << "Ending test function createPartitionManagerForSubGroup" << endl;

   // cleaning
   splitor.clearInternal();
   m_pm.clear();
   AOMD::AOMD_Util::Instance()->deleteMeshDataId(tag_level);

   // wait every destruction done
   MPI_Barrier(MPI_COMM_WORLD);

   return xMPIEnv::finalize();
}
