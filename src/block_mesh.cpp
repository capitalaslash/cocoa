#include "problem/problem_oforg.hpp"

// stl
#include <memory>

// openfoam.org - blockMesh
#include <OFstream.H>
#include <attachPolyTopoChanger.H>
#include <blockMesh.H>
#include <polyTopoChange.H>
#include <slidingInterface.H>
#include <systemDict.H>

int runBlockMesh(Foam::argList & args)
{
  // TODO: options should be defined before the creation of the args objest
  Foam::argList::noParallel();

  Foam::argList::addOption(
      "dict", "file", "read control dictionary from specified location");

  Foam::argList::addBoolOption(
      "blockTopology", "write block edges and centres as .obj files");
  Foam::argList::addBoolOption("noClean", "keep the existing files in the polyMesh");

  Foam::argList::addNote(
      "Block description\n"
      "\n"
      "  For a given block, the correspondence between the ordering of\n"
      "  vertex labels and face labels is shown below.\n"
      "  For vertex numbering in the sequence 0 to 7 (block, centre):\n"
      "    faces 0 (f0) and 1 are left and right, respectively;\n"
      "    faces 2 and 3 are front and back; \n"
      "    and faces 4 and 5 are bottom and top::\n"
      "\n"
      "                 7 ---- 6\n"
      "            f5   |\\     |\\   f3\n"
      "            |    | 4 ---- 5   \\\n"
      "            |    3 |--- 2 |    \\\n"
      "            |     \\|     \\|    f2\n"
      "            f4     0 ---- 1\n"
      "\n"
      "       Z         f0 ----- f1\n"
      "       |  Y\n"
      "       | /\n"
      "       O --- X\n");

  if (!args.checkRootCase())
  {
    Foam::FatalError.exit();
  }

  Foam::Info << "Create time\n" << Foam::endl;

  Foam::Time runTime(Foam::Time::controlDictName, args);

  const Foam::word dictName("blockMeshDict");

  Foam::word regionName;
  Foam::word regionPath;

  // Check if the region is specified otherwise mesh the default region
  if (args.optionReadIfPresent("region", regionName, Foam::polyMesh::defaultRegion))
  {
    Foam::Info << Foam::nl << "Generating mesh for region " << regionName << Foam::endl;
    regionPath = regionName;
  }

  if (!args.optionFound("noClean"))
  {
    Foam::fileName polyMeshPath(
        runTime.path() / runTime.constant() / regionPath / Foam::polyMesh::meshSubDir);

    if (exists(polyMeshPath))
    {
      if (exists(polyMeshPath / dictName))
      {
        Foam::Info << "Not deleting polyMesh directory " << Foam::nl << "    "
                   << polyMeshPath << Foam::nl << "    because it contains " << dictName
                   << Foam::endl;
      }
      else
      {
        Foam::Info << "Deleting polyMesh directory" << Foam::nl << "    "
                   << polyMeshPath << Foam::endl;
        Foam::rmDir(polyMeshPath);
      }
    }
  }

  Foam::typeIOobject<Foam::IOdictionary> meshDictIO(
      Foam::systemDictIO(dictName, args, runTime, regionName));

  if (!meshDictIO.headerOk())
  {
    FatalErrorInFunction << "Cannot find file " << meshDictIO.relativeObjectPath()
                         << Foam::nl << Foam::exit(Foam::FatalError);
  }

  Foam::Info << "Creating block mesh from\n    " << meshDictIO.relativeObjectPath()
             << Foam::endl;

  Foam::IOdictionary meshDict(meshDictIO);
  Foam::blockMesh blocks(meshDict, regionName);

  if (args.optionFound("blockTopology"))
  {
    // Write mesh as edges.
    {
      Foam::fileName objMeshFile("blockTopology.obj");

      Foam::OFstream str(runTime.path() / objMeshFile);

      Foam::Info << Foam::nl << "Dumping block structure as Lightwave obj format"
                 << " to " << objMeshFile << Foam::endl;

      blocks.writeTopology(str);
    }

    // Write centres of blocks
    {
      Foam::fileName objCcFile("blockCentres.obj");

      Foam::OFstream str(runTime.path() / objCcFile);

      Foam::Info << Foam::nl << "Dumping block centres as Lightwave obj format"
                 << " to " << objCcFile << Foam::endl;

      const Foam::polyMesh & topo = blocks.topology();

      const Foam::pointField & cellCentres = topo.cellCentres();

      forAll(cellCentres, celli)
      {
        // point cc = b.blockShape().centre(b.points());
        const Foam::point & cc = cellCentres[celli];

        str << "v " << cc.x() << ' ' << cc.y() << ' ' << cc.z() << Foam::nl;
      }
    }

    Foam::Info << Foam::nl << "end" << Foam::endl;

    return 0;
  }

  Foam::Info << Foam::nl << "Creating polyMesh from blockMesh" << Foam::endl;

  Foam::word defaultFacesName = "defaultFaces";
  Foam::word defaultFacesType = Foam::emptyPolyPatch::typeName;
  Foam::polyMesh mesh(
      Foam::IOobject(regionName, runTime.constant(), runTime),
      clone(blocks.points()), // could we re-use space?
      blocks.cells(),
      blocks.patches(),
      blocks.patchNames(),
      blocks.patchDicts(),
      defaultFacesName,
      defaultFacesType);

  // Read in a list of dictionaries for the merge patch pairs
  if (meshDict.found("mergePatchPairs"))
  {
    Foam::List<Foam::Pair<Foam::word>> mergePatchPairs(
        meshDict.lookup("mergePatchPairs"));

    if (mergePatchPairs.size())
    {
      Foam::Info << "Creating merge patch pairs" << Foam::nl << Foam::endl;

      // Create and add point and face zones and mesh modifiers
      Foam::List<std::unique_ptr<Foam::pointZone>> pz(mergePatchPairs.size());
      Foam::List<std::unique_ptr<Foam::faceZone>> fz(3 * mergePatchPairs.size());
      Foam::List<Foam::cellZone *> cz(0);

      forAll(mergePatchPairs, pairI)
      {
        const Foam::word mergeName(
            mergePatchPairs[pairI].first() + mergePatchPairs[pairI].second() +
            Foam::name(pairI));

        pz[pairI].reset(new Foam::pointZone(
            mergeName + "CutPointZone", Foam::labelList(0), 0, mesh.pointZones()));

        // Master patch
        const Foam::word masterPatchName(mergePatchPairs[pairI].first());
        const Foam::polyPatch & masterPatch = mesh.boundaryMesh()[masterPatchName];

        Foam::labelList isf(masterPatch.size());

        forAll(isf, i) { isf[i] = masterPatch.start() + i; }

        fz[3 * pairI].reset(new Foam::faceZone(
            mergeName + "MasterZone",
            isf,
            Foam::boolList(masterPatch.size(), false),
            0,
            mesh.faceZones()));

        // Slave patch
        const Foam::word slavePatchName(mergePatchPairs[pairI].second());
        const Foam::polyPatch & slavePatch = mesh.boundaryMesh()[slavePatchName];

        Foam::labelList osf(slavePatch.size());

        forAll(osf, i) { osf[i] = slavePatch.start() + i; }

        fz[3 * pairI + 1].reset(new Foam::faceZone(
            mergeName + "SlaveZone",
            osf,
            Foam::boolList(slavePatch.size(), false),
            1,
            mesh.faceZones()));

        // Add empty zone for cut faces
        fz[3 * pairI + 2].reset(new Foam::faceZone(
            mergeName + "CutFaceZone",
            Foam::labelList(0),
            Foam::boolList(0, false),
            2,
            mesh.faceZones()));
      } // end of all merge pairs

      Foam::Info << "Adding point and face zones" << Foam::endl;
      Foam::List<Foam::pointZone *> pzPtrs(mergePatchPairs.size());
      Foam::List<Foam::faceZone *> fzPtrs(3 * mergePatchPairs.size());
      for (auto k = 0; k < mergePatchPairs.size(); k++)
      {
        pzPtrs[k] = pz[k].get();
        fzPtrs[3 * k + 0] = fz[3 * k + 0].get();
        fzPtrs[3 * k + 1] = fz[3 * k + 1].get();
        fzPtrs[3 * k + 2] = fz[3 * k + 2].get();
      }
      mesh.addZones(pzPtrs, fzPtrs, cz);

      Foam::Info << "Creating attachPolyTopoChanger" << Foam::endl;
      Foam::attachPolyTopoChanger polyMeshAttacher(mesh);
      polyMeshAttacher.setSize(mergePatchPairs.size());

      forAll(mergePatchPairs, pairI)
      {
        const Foam::word mergeName(
            mergePatchPairs[pairI].first() + mergePatchPairs[pairI].second() +
            Foam::name(pairI));

        // Add the sliding interface mesh modifier
        auto slidingInterface =
            std::unique_ptr<Foam::slidingInterface>{new Foam::slidingInterface{
                "couple" + Foam::name(pairI),
                pairI,
                polyMeshAttacher,
                mergeName + "MasterZone",
                mergeName + "SlaveZone",
                mergeName + "CutPointZone",
                mergeName + "CutFaceZone",
                mergePatchPairs[pairI].first(),
                mergePatchPairs[pairI].second(),
                Foam::slidingInterface::INTEGRAL, // always integral
                false,
                Foam::intersection::algorithm::visible}};
        polyMeshAttacher.set(pairI, slidingInterface.get());
      }

      polyMeshAttacher.attach(true);
    }
  }
  else
  {
    Foam::Info << Foam::nl << "There are no merge patch pairs edges" << Foam::endl;
  }

  // Set any cellZones (note: cell labelling unaffected by above
  // mergePatchPairs)

  Foam::label nZones = blocks.numZonedBlocks();

  if (nZones > 0)
  {
    Foam::Info << Foam::nl << "Adding cell zones" << Foam::endl;

    // Map from zoneName to cellZone index
    Foam::HashTable<Foam::label> zoneMap(nZones);

    // Cells per zone.
    Foam::List<Foam::DynamicList<Foam::label>> zoneCells(nZones);

    // Running cell counter
    Foam::label celli = 0;

    // Largest zone so far
    Foam::label freeZoneI = 0;

    forAll(blocks, blockI)
    {
      const Foam::block & b = blocks[blockI];
      const Foam::List<Foam::FixedList<Foam::label, 8>> blockCells = b.cells();
      const Foam::word & zoneName = b.zoneName();

      if (zoneName.size())
      {
        Foam::HashTable<Foam::label>::const_iterator iter = zoneMap.find(zoneName);

        Foam::label zoneI;

        if (iter == zoneMap.end())
        {
          zoneI = freeZoneI++;

          Foam::Info << "    " << zoneI << '\t' << zoneName << Foam::endl;

          zoneMap.insert(zoneName, zoneI);
        }
        else
        {
          zoneI = iter();
        }

        forAll(blockCells, i) { zoneCells[zoneI].append(celli++); }
      }
      else
      {
        celli += blockCells.size();
      }
    }

    Foam::List<std::unique_ptr<Foam::cellZone>> cz(zoneMap.size());

    forAllConstIter(Foam::HashTable<Foam::label>, zoneMap, iter)
    {
      Foam::label zoneI = iter();

      cz[zoneI].reset(new Foam::cellZone(
          iter.key(), zoneCells[zoneI].shrink(), zoneI, mesh.cellZones()));
    }

    mesh.pointZones().setSize(0);
    mesh.faceZones().setSize(0);
    mesh.cellZones().setSize(0);
    Foam::List<Foam::cellZone *> czPtrs(zoneMap.size());
    for (int k = 0; k < zoneMap.size(); k++)
      czPtrs[k] = cz[k].get();
    mesh.addZones(
        Foam::List<Foam::pointZone *>(0), Foam::List<Foam::faceZone *>(0), czPtrs);
  }

  // Detect any cyclic patches and force re-ordering of the faces
  {
    const Foam::polyPatchList & patches = mesh.boundaryMesh();
    bool hasCyclic = false;
    forAll(patches, patchi)
    {
      if (Foam::isA<Foam::cyclicPolyPatch>(patches[patchi]))
      {
        hasCyclic = true;
        break;
      }
    }

    if (hasCyclic)
    {
      Foam::Info << Foam::nl << "Detected cyclic patches; ordering boundary faces"
                 << Foam::endl;
      const Foam::word oldInstance = mesh.instance();
      Foam::polyTopoChange meshMod(mesh);
      meshMod.changeMesh(mesh, false);
      mesh.setInstance(oldInstance);
    }
  }

  // Set the precision of the points data to 10
  Foam::IOstream::defaultPrecision(Foam::max(10u, Foam::IOstream::defaultPrecision()));

  Foam::Info << Foam::nl << "Writing polyMesh" << Foam::endl;
  mesh.removeFiles();
  if (!mesh.write())
  {
    FatalErrorInFunction << "Failed writing polyMesh." << Foam::exit(Foam::FatalError);
  }

  // Write summary
  {
    const Foam::polyPatchList & patches = mesh.boundaryMesh();

    Foam::Info << "----------------" << Foam::nl << "Mesh Information" << Foam::nl
               << "----------------" << Foam::nl << "  "
               << "boundingBox: " << Foam::boundBox(mesh.points()) << Foam::nl << "  "
               << "nPoints: " << mesh.nPoints() << Foam::nl << "  "
               << "nCells: " << mesh.nCells() << Foam::nl << "  "
               << "nFaces: " << mesh.nFaces() << Foam::nl << "  "
               << "nInternalFaces: " << mesh.nInternalFaces() << Foam::nl;

    Foam::Info << "----------------" << Foam::nl << "Patches" << Foam::nl
               << "----------------" << Foam::nl;

    forAll(patches, patchi)
    {
      const Foam::polyPatch & p = patches[patchi];

      Foam::Info << "  " << "patch " << patchi << " (start: " << p.start()
                 << " size: " << p.size() << ") name: " << p.name() << Foam::nl;
    }
  }

  Foam::Info << "\nEnd blockMesh\n" << Foam::endl;

  return 0;
}
