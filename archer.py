import os, collections

from SBI                       import SBIglobals
from SBI.src.structure         import PDB
from SBI.src.beans             import Path, File
from SBI.src.structure.protein import Arch

def build_archs(sourcepdb, 
                chain_select        = None, 
                limit_internal_ss   = 100, 
                limit_distance      = False, 
                allowed_gaps        = 0, 
                clean_dssp_files    = True,
                workdir             = os.getcwd()):

    newpdb = PDB(sourcepdb, dehydrate = True)

    if chain_select is None: chain_select = newpdb.chain_identifiers
    else: 
        if not isinstance(chain_select, collections.Iterable): chain_select = set([chain_select])

    for working_chain in newpdb.proteins:
        if working_chain.chain not in chain_select: continue

        working_chain.calculate_dssp(os.path.join(workdir, 'tmppdb.' + working_chain.globalID), 
                                     os.path.join(workdir, 'tmpdssp.' + working_chain.globalID),
                                     clean_dssp_files)

        working_chain.calculate_archs(limit_internal_ss = limit_internal_ss, 
                                      limit_distance    = limit_distance, 
                                      allowed_gaps      = allowed_gaps)

    return newpdb

def sortarchs(inputdir, outputdir):
    
    archsdir              = outputdir
    Path.mkdir(archsdir)
    sorted_archs          = {}
    loop_file_name        = os.path.join(archsdir, 'ArchDB.{0}.db')
    loop_split_file_name  = os.path.join(archsdir, 'ArchDB.{0}.{1:02d}-{2:02d}.db')
    sections_ini          = [ 0, 4, 7,14,21]
    sections_end          = [ 4, 6,13,20, 0]
    for archfile in Path.list_files(root = inputdir, pattern = '*.archObj'):
        filename = os.path.basename(archfile)
        data     = filename.split('_')
        length   = int(data[0])
        archtype = data[1] 
        sorted_archs.setdefault(archtype,{}).setdefault(length,[])
        sorted_archs[archtype][length].append(archfile)
    
    for archtype in sorted_archs:
        SBIglobals.alert('verbose', None, "ARCHS: " + archtype + "\n")
        fd  = File(loop_file_name.format(archtype), 'w')
        fdp = []
        for x in range(len(sections_ini)):
            fdp.append(File(loop_split_file_name.format(archtype, sections_ini[x], sections_end[x]), 'w'))
        for length in sorted(sorted_archs[archtype]):
            SBIglobals.alert('verbose', None, '\t{0}'.format(length))
            for archfile in sorted_archs[archtype][length]:
                SBIglobals.alert('verbose', None, '\t\t{0}\n'.format(archfile))
                nsp = Arch.load(archfile)
                fd.descriptor.write(nsp.archtype_format() + "\n")
                for x in range(len(fdp)):
                    if length >= sections_ini[x] and (sections_end[x] == 0 or length <= sections_end[x]):
                        fdp[x].descriptor.write(nsp.archtype_format() + "\n")
        fd.close()
        for x in range(len(fdp)):
            fdp[x].close()

def join_PDB_archs(archlist, prefix, js=False):

    chain_id_number = 65 #chr(65) = A

    pdbout    = open(prefix + '.pdb', 'w')
    if js: 
        jsout = open(prefix + '.js',  'w')
        jsout.write('var {0}="'.format(prefix.split('/')[-1]))

    chains = []
    for pdb_file in archlist:
        newpdb = PDB(pdb_file, header = True)
        for ss in newpdb.header.secondary_structures:
            ss.change_chain(chr(chain_id_number))
            pdbout.write(str(ss) + "\n")
            if js: jsout.write(str(ss) + "\\n")
        chains.append(newpdb.chains[0])
        chains[-1].chain = chr(chain_id_number)
        chain_id_number += 1

    for chain in chains:
        pdbout.write(chain.PDB_format() + "\n")
        if js: jsout.write(chain.PDB_format().replace('\n','\\n') + "\\n")

    pdbout.close()
    if js:
        jsout.write('";') 
        jsout.close()