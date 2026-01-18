OBJ=obj
BIN=bin
OUT=output
SDF=sdf
TMP=tmp

DIRS=$(OBJ) $(BIN) $(OUT) $(SDF) $(TMP)
OBJS=$(OBJ)/misc.o $(OBJ)/point.o $(OBJ)/atom.o $(OBJ)/intera.o $(OBJ)/molecule.o $(OBJ)/aminoacid.o \
	$(OBJ)/protein.o $(OBJ)/moiety.o $(OBJ)/conj.o $(OBJ)/progress.o
DOBJ=$(OBJ)/dynamic.o $(OBJ)/reshape.o $(OBJ)/scoring.o $(OBJ)/search.o $(OBJ)/cavity.o $(OBJ)/soft.o $(OBJ)/appear.o
TSTS=test/point_test test/atom_test test/molecule_test test/pi_stack_test test/mol_assem_test test/aniso_test test/amino_test \
	  test/protein_test test/backbone_test test/bond_rotation_test test/moiety_test test/ameliorate_test \
	  test/flexion_test test/histidine_test test/ring_test test/eclipsing_test test/mcoord_test test/vdw_vertex_test \
	  test/ageo_test test/chirality_test test/bb_test test/solvent_test test/multimer_test test/inte_test test/conj_test
APPS=$(BIN)/aromadock $(BIN)/phew $(BIN)/ic $(BIN)/qc $(BIN)/protseq $(BIN)/molsurf $(BIN)/olfactophore \
	 $(BIN)/scorpion $(BIN)/ramachandran $(BIN)/ringflip $(BIN)/cavity_search $(BIN)/cavity_fit
all: $(DIRS) \
	 $(OBJS) \
	 $(DOBJ) \
	 $(TSTS) \
	 $(APPS)
code: $(DIRS) $(OBJS) $(DOBJ) $(TSTS) $(APPS)
apps: $(APPS)
aromadock: $(DIRS) $(OBJS) $(DOBJ) $(BIN)/aromadock
phew: $(DIRS) $(OBJS) $(DOBJ) $(BIN)/phew
ic: $(DIRS) $(OBJS) $(DOBJ) $(BIN)/ic

CPL=g++

# Common flags for all modes
CFLAGS=-ffast-math -Wwrite-strings -fextended-identifiers -std=c++14

# Default CFLAGS - release mode
CFLAGS+=-O3

# Debug CFLAGS - allows gdb, valgrind
# CFLAGS+=-g

# For gprof
# example command line:
# gprof bin/aromadock gmon.out > aromadock.output
# CFLAGS+=-g -pg

# For code coverage instrumentation, switch to these CFLAGS (slower performance):
# CFLAGS+=-g -fprofile-arcs -ftest-coverage

clean:
	rm -f $(OBJ)/*.o $(BIN)/*

$(OBJ):
	if [ ! -f $(OBJ) ]; then mkdir -p $(OBJ); fi

$(BIN):
	if [ ! -f $(BIN) ]; then mkdir -p $(BIN); fi

$(OUT):
	if [ ! -f $(OUT) ]; then mkdir -p $(OUT); fi

$(SDF):
	if [ ! -f $(SDF) ]; then mkdir -p $(SDF); fi

$(TMP):
	if [ ! -f $(TMP) ]; then mkdir -p $(TMP); fi


# Classes

$(OBJ)/misc.o: src/classes/misc.h src/classes/misc.cpp src/classes/constants.h makefile
	$(CPL) -c src/classes/misc.cpp -o $(OBJ)/misc.o $(CFLAGS)

$(OBJ)/point.o: src/classes/point.h src/classes/point.cpp $(OBJ)/misc.o src/classes/constants.h
	$(CPL) -c src/classes/point.cpp -o $(OBJ)/point.o $(CFLAGS)

$(OBJ)/atom.o: src/classes/atom.h src/classes/atom.cpp $(OBJ)/point.o
	$(CPL) -c src/classes/atom.cpp -o $(OBJ)/atom.o $(CFLAGS)

$(OBJ)/conj.o: src/classes/conj.h src/classes/conj.cpp $(OBJ)/atom.o
	$(CPL) -c src/classes/conj.cpp -o $(OBJ)/conj.o $(CFLAGS)

$(OBJ)/intera.o: src/classes/intera.h src/classes/intera.cpp $(OBJ)/conj.o
	$(CPL) -c src/classes/intera.cpp -o $(OBJ)/intera.o $(CFLAGS)

$(OBJ)/molecule.o: src/classes/molecule.h src/classes/molecule.cpp $(OBJ)/intera.o
	$(CPL) -c src/classes/molecule.cpp -o $(OBJ)/molecule.o $(CFLAGS)

$(OBJ)/aminoacid.o: src/classes/aminoacid.h src/classes/aminoacid.cpp $(OBJ)/molecule.o
	$(CPL) -c src/classes/aminoacid.cpp -o $(OBJ)/aminoacid.o $(CFLAGS)

$(OBJ)/protein.o: src/classes/protein.h src/classes/protein.cpp $(OBJ)/aminoacid.o
	$(CPL) -c src/classes/protein.cpp -o $(OBJ)/protein.o $(CFLAGS)

$(OBJ)/reshape.o: src/classes/reshape.h src/classes/reshape.cpp $(OBJ)/protein.o
	$(CPL) -c src/classes/reshape.cpp -o $(OBJ)/reshape.o $(CFLAGS)

$(OBJ)/search.o: src/classes/search.h src/classes/search.cpp $(OBJ)/protein.o
	$(CPL) -c src/classes/search.cpp -o $(OBJ)/search.o $(CFLAGS)

$(OBJ)/cavity.o: src/classes/cavity.h src/classes/cavity.cpp $(OBJ)/protein.o
	$(CPL) -c src/classes/cavity.cpp -o $(OBJ)/cavity.o $(CFLAGS)

$(OBJ)/soft.o: src/classes/soft.h src/classes/soft.cpp $(OBJ)/protein.o
	$(CPL) -c src/classes/soft.cpp -o $(OBJ)/soft.o $(CFLAGS)

$(OBJ)/appear.o: src/classes/appear.h src/classes/appear.cpp $(OBJ)/protein.o
	$(CPL) -c src/classes/appear.cpp -o $(OBJ)/appear.o $(CFLAGS)

$(OBJ)/dynamic.o: src/classes/dynamic.h src/classes/dynamic.cpp $(OBJ)/protein.o
	$(CPL) -c src/classes/dynamic.cpp -o $(OBJ)/dynamic.o $(CFLAGS)

$(OBJ)/moiety.o: src/classes/moiety.h src/classes/moiety.cpp $(OBJ)/molecule.o
	$(CPL) -c src/classes/moiety.cpp -o $(OBJ)/moiety.o $(CFLAGS)

$(OBJ)/scoring.o: src/classes/scoring.h src/classes/scoring.cpp $(OBJ)/search.o
	$(CPL) -c src/classes/scoring.cpp -o $(OBJ)/scoring.o $(CFLAGS)

$(OBJ)/progress.o: src/classes/progress.h src/classes/progress.cpp $(OBJ)/misc.o
	$(CPL) -c src/classes/progress.cpp -o $(OBJ)/progress.o $(CFLAGS)


# Tests

test/point_test: src/test/point_test.cpp $(OBJ)/point.o
	$(CPL) src/test/point_test.cpp $(OBJ)/point.o $(OBJ)/misc.o -o test/point_test $(CFLAGS)

test/atom_test: src/test/atom_test.cpp $(OBJ)/point.o $(OBJ)/atom.o
	$(CPL) src/test/atom_test.cpp $(OBJ)/misc.o $(OBJ)/atom.o $(OBJ)/point.o -o test/atom_test $(CFLAGS)

test/molecule_test: src/test/molecule_test.cpp $(OBJS)
	$(CPL) src/test/molecule_test.cpp $(OBJS) -o test/molecule_test $(CFLAGS)

test/ring_test: src/test/ring_test.cpp $(OBJ)/molecule.o
	$(CPL) src/test/ring_test.cpp $(OBJS) -o test/ring_test $(CFLAGS)

test/pi_stack_test: src/test/pi_stack_test.cpp $(OBJS)
	$(CPL) src/test/pi_stack_test.cpp $(OBJS) -o test/pi_stack_test $(CFLAGS)

test/aniso_test: src/test/aniso_test.cpp $(OBJS)
	$(CPL) src/test/aniso_test.cpp $(OBJS) -o test/aniso_test $(CFLAGS)

test/inte_test: src/test/inte_test.cpp $(OBJS)
	$(CPL) src/test/inte_test.cpp $(OBJS) -o test/inte_test $(CFLAGS)

test/conj_test: src/test/conj_test.cpp $(OBJS)
	$(CPL) src/test/conj_test.cpp $(OBJS) -o test/conj_test $(CFLAGS)

test/moiety_test: src/test/moiety_test.cpp $(OBJS)
	$(CPL) src/test/moiety_test.cpp $(OBJS) -o test/moiety_test $(CFLAGS)

test/ameliorate_test: src/test/ameliorate_test.cpp $(OBJS)
	$(CPL) src/test/ameliorate_test.cpp $(OBJS) -o test/ameliorate_test $(CFLAGS)

test/mol_assem_test: src/test/mol_assem_test.cpp $(OBJS)
	$(CPL) src/test/mol_assem_test.cpp $(OBJS) -o test/mol_assem_test $(CFLAGS)

test/eclipsing_test: src/test/eclipsing_test.cpp $(OBJS)
	$(CPL) src/test/eclipsing_test.cpp $(OBJS) -o test/eclipsing_test $(CFLAGS)

test/amino_test: src/test/amino_test.cpp $(OBJS) $(OBJ)/aminoacid.o
	$(CPL) src/test/amino_test.cpp $(OBJS) -o test/amino_test $(CFLAGS)

test/protein_test: src/test/protein_test.cpp $(OBJS) $(OBJ)/aminoacid.o $(OBJ)/protein.o
	$(CPL) src/test/protein_test.cpp $(OBJS) -o test/protein_test $(CFLAGS)

test/ageo_test: src/test/ageo_test.cpp $(OBJS) $(OBJ)/aminoacid.o $(OBJ)/protein.o
	$(CPL) src/test/ageo_test.cpp $(OBJS) -o test/ageo_test $(CFLAGS)

test/chirality_test: src/test/chirality_test.cpp $(OBJ)/molecule.o
	$(CPL) src/test/chirality_test.cpp $(OBJS) -o test/chirality_test $(CFLAGS)

test/backbone_test: src/test/backbone_test.cpp $(OBJS) $(OBJ)/aminoacid.o $(OBJ)/protein.o
	$(CPL) src/test/backbone_test.cpp $(OBJS) -o test/backbone_test $(CFLAGS)

test/bond_rotation_test: src/test/bond_rotation_test.cpp $(OBJS) $(OBJ)/molecule.o
	$(CPL) src/test/bond_rotation_test.cpp $(OBJS) -o test/bond_rotation_test $(CFLAGS)

test/flexion_test: src/test/flexion_test.cpp $(OBJS)
	$(CPL) src/test/flexion_test.cpp $(OBJS) -o test/flexion_test $(CFLAGS)

test/histidine_test: src/test/histidine_test.cpp $(OBJS)
	$(CPL) src/test/histidine_test.cpp $(OBJS) -o test/histidine_test $(CFLAGS)

test/bb_test: src/test/bb_test.cpp $(OBJS) $(DOBJ)
	$(CPL) src/test/bb_test.cpp $(OBJS) $(DOBJ) -o test/bb_test $(CFLAGS)

test/solvent_test: src/test/solvent_test.cpp $(OBJS)
	$(CPL) src/test/solvent_test.cpp $(OBJS) -o test/solvent_test $(CFLAGS)

test/multimer_test: src/test/multimer_test.cpp $(OBJS)
	$(CPL) src/test/multimer_test.cpp $(OBJS) -o test/multimer_test $(CFLAGS)

test/mcoord_test: src/test/mcoord_test.cpp $(OBJS)
	$(CPL) src/test/mcoord_test.cpp $(OBJS) -o test/mcoord_test $(CFLAGS)

test/vdw_vertex_test: src/test/vdw_vertex_test.cpp $(OBJS)
	$(CPL) src/test/vdw_vertex_test.cpp $(OBJS) -o test/vdw_vertex_test $(CFLAGS)

test/eclipse: src/test/eclipse.cpp $(OBJS) $(DOBJ)
	$(CPL) src/test/eclipse.cpp $(OBJS) -o test/eclipse $(CFLAGS)


# Apps

$(BIN)/aromadock: src/aromadock.cpp $(OBJS) $(DOBJ)
	$(CPL) src/aromadock.cpp $(OBJS) $(DOBJ) -o $(BIN)/aromadock $(CFLAGS)

$(BIN)/phew: src/phew.cpp $(OBJS) $(DOBJ)
	$(CPL) src/phew.cpp $(OBJS) $(DOBJ) -o $(BIN)/phew $(CFLAGS)

$(BIN)/cavity_search: src/cavity_search.cpp $(OBJS) $(DOBJ)
	$(CPL) src/cavity_search.cpp $(OBJS) $(DOBJ) -o $(BIN)/cavity_search $(CFLAGS)

$(BIN)/cavity_fit: src/cavity_fit.cpp $(OBJS) $(DOBJ)
	$(CPL) src/cavity_fit.cpp $(OBJS) $(DOBJ) -o $(BIN)/cavity_fit $(CFLAGS)

$(BIN)/protseq: src/protseq.cpp $(OBJS) $(DOBJ)
	$(CPL) src/protseq.cpp $(OBJS) $(DOBJ) -o $(BIN)/protseq $(CFLAGS)

$(BIN)/ic: src/ic.cpp $(OBJS) $(DOBJ)
	$(CPL) src/ic.cpp $(OBJS) $(DOBJ) -o $(BIN)/ic $(CFLAGS)

$(BIN)/qc: src/qc.cpp $(OBJS) $(DOBJ)
	$(CPL) src/qc.cpp $(OBJS) $(DOBJ) -o $(BIN)/qc $(CFLAGS)

$(BIN)/scorpion: src/scorpion.cpp $(OBJS) $(DOBJ)
	$(CPL) src/scorpion.cpp $(OBJS) $(DOBJ) -o $(BIN)/scorpion $(CFLAGS)

$(BIN)/olfactophore: src/olfactophore.cpp $(OBJS) $(DOBJ)
	$(CPL) src/olfactophore.cpp $(OBJS) $(DOBJ) -o $(BIN)/olfactophore $(CFLAGS)

$(BIN)/ramachandran: src/ramachandran.cpp $(OBJS) $(DOBJ)
	$(CPL) src/ramachandran.cpp $(OBJS) $(DOBJ) -o $(BIN)/ramachandran $(CFLAGS)

$(BIN)/ringflip: src/ringflip.cpp $(OBJS)
	$(CPL) src/ringflip.cpp $(OBJS) -o $(BIN)/ringflip $(CFLAGS)

$(BIN)/molsurf: src/molsurf.cpp $(OBJS)
	$(CPL) src/molsurf.cpp $(OBJS) -o $(BIN)/molsurf $(CFLAGS)
