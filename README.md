
https://www.primaryodors.org

# Dependencies

g++ make openbabel python3 python3-pybel python3-biopython python3-natsort librdkit1 python3-rdkit rdkit-data php php-curl php-gd


# Performing Docks

First navigate to the AROMA root folder and run `make`. Note all command line tools must be run from this project root folder.

Next, type `./dock.sh` with the receptor as the first argument and the odorant as the second argument. This will attempt to dock the ligand
inside the active and inactive models of the receptor.

In lieu of a receptor, you may specify:
- `all` meaning to run for all receptors;
- `emp` meaning all empirical pairs for the odorant(s);
- `ago` meaning all receptors where each odorant is a known agonist;
- `e1` ~ `e100` meaning only receptors whose rate of expression in the human population is known to be at least the number (percentage) specified.

Regarding expression rates, it was found by Verbeurgt et al (2014) that each olfactory receptor is expressed in a subset of individuals.
They searched whole olfactory mucosa of 26 cadavers for RNA signatures of the ORs and found that some receptors were expressed in all subjects,
while others ranged from common to rare in their expression rates. Therefore, AromaDock allows filtering for only those receptors expressed
at or above some user specified percentage. Interestingly, no subject in the study expressed the full set of all ~400 receptors.

In lieu of an odorant, you may specify:
- `all` meaning all odorants in the database;
- `emp` meaning all empirical pairs for the receptor(s);
- `ago` meaning all known agonists of the receptor(s);
- a perceptual note from the odorants in the database, such as `citrus`, `dairy`, or `acrid`.


# Web Application

You may optionally host your own AROMA web interface.

To enable the web app:
- Either set up a local web server or checkout this repository in a folder on a web host.
- Make sure your server has the `php`, `php-curl`, `php-gd`, and `openbabel` packages installed.
- After installing `php-curl`, it's important to restart the web service e.g. `sudo apache2ctl -k restart`.
- Then open the `www/symlink.sh` file in a text editor, make sure the destination folder is correct (by default it will show `/var/www/html/`
  which is usually correct for Apache2 installations), make sure you have write permissions in the 
  folder (or use `sudo`), and execute `www/symlink.sh` i n a command line.
- The `data` and `www/assets` folders and all contents must also be recursively made writable by the web user.
- If on a local server, you will now have an instance of the web app at http://127.0.0.1/aroma/ whereas if you are using a web host
  then you may have to configure your hosting to point one of your registered domains or subfolders to the `aroma/www` folder.

If you get a 403 Forbidden error, please make sure that every containing folder of the `aroma/www` folder has public execute access.


# Homology Modeling

All of the active models under the `pdbs/` folder have been created with the `hm/dohm.php` script. Many have also undergone post-processing steps
with other utilities such as `bin/ic`. To create your own active-state homology models, simply run `php -f hm/dohm.php OR###` from the `aroma/`
project root folder, replacing `OR###` with the ID of the receptor you wish to model. WARNING: This action WILL OVERWRITE existing active models.

You can also adjust an active-state model if a known agonist fails to dock, by running `hm/fixfail.py` from the `aroma/` folder. The receptor ID 
and the docked PDB model are the only required parameters. WARNING: This action WILL OVERWRITE existing active models. An example usage might be:

`python3 hm/fixfail.py OR1A2 out/OR1/OR1A2/OR1A2~nerolidol.active.model1.pdb`

This example will recreate the active-state OR1A2 model using the output PDB as the template and the docked ligand from the PDB as a rigid body.
See https://salilab.org/modeller/9v1/manual/node94.html for more information on rigid bodies. After running `hm/fixfail.py`, it is a good idea to
retry the dock that previously failed in order to check whether the new model is satisfactory. If not, `hm/fixfail.py` can be rerun an unlimited
number of times, retrying the dock after each time.


# Tools

Note: All AROMA tools are designed to be run in a command prompt from the `aroma/` project root folder. Running them from any other working
directory may result in errors, file-not-found messages, exceptions, and undefined behavior.

`bin/aromadock`
The molecular docker.

`bin/cavity_fit`
A tool for placing odorant molecules inside cavities identified by bin/cavity_search.

`bin/cavity_search`
A tool for scanning proteins for ligand binding cavities.

`bin/ic`
A tool for scanning proteins for internal contacts (e.g. inter-residue hydrogen bonds, salt bridges, etc.) as well as limited protein editing,
internal clash minimization, and hydrogen position optimization.

`bin/molsurf`
A tool for estimating the surface are of molecules.

`bin/olfactophore`
A tool for superimposing odorant molecules in space in order to identify the spatial and electrostatic features of collections of odorants,
such as those with a common aroma note or the agonists of a given receptor.

`bin/phew`
Protein Editing and Hydrogenation Workspace. A script interpreter for more advanced protein editing.

`bin/protseq`
Reads a PDB file and outputs the amino acid sequence.

`bin/qc`
Reads a multi-strand PDB file and scans for contacts (e.g. hydrogen bonds, salt bridges, etc.) between two indicated strands.

`bin/ramachandran`
Outputs Ramachandran diagrams of PDB files.

`bin/ringflip`
Edits SDF files, performing "flips" of ring atoms, for example to toggle between chair and boat forms of cyclohexane or hexoses, etc.

`bin/scorpion`
Formerly ScorePDB, a tool that reads a ligand-bound PDB and re-runs the same scoring code as AromaDock to generate .dock format output.

`data/btree.php`
A script to rebuild the phylogenic tree of receptors. WARNING: This script WILL OVERWRITE data/receptor.json and data/treenodes.json!

`data/comparison.py`
A tool for comparing the amino acids at one or more Ballesteros-Weinstein residues across multiple receptors.

`data/conserved.py`
A tool for finding the distribution of aminos at a single Ballesteros-Weinstein position across all receptors.

`data/dyncenter.py`
Used by the docking script. Finds the correct binding site residues (BSRs) for a given receptor+ligand pair.

`data/generate_fasta.php`
Outputs FASTA format for one receptor or all receptors.

`data/gensdf.py`
A tool for adding new odorants to the database or generating 3D SDF structures of existng odorants.

`data/measurebw.py`
A tool for measuring distances between residues in a PDB.

`data/odordupes.php`
A tool for finding possible duplicate entries in the odorants database.

`data/sequence_update.php`
Updates the receptors database with information from the data/sequences_aligned.txt file. WARNING: This script WILL OVERWRITE data/receptor.json.

`hm/dohm.php`
The homology modeling utility. Requires MODELLER in order to run. WARNING: This script WILL OVERWRITE receptor PDB models.

`hm/fixfail.py`
A script for recreating new active-state homology models based off of failed agonist docks. WARNING: This script WILL OVERWRITE existing models.

`dock.sh`
A shell wrapper for run.py.

`reec.py`
A work in progress. Receptor Energetic Efficacy Calculator. Ultimately a script to examine dock results and make predictions about receptor responses
to odorants.

`run.py`
Accepts a receptor filter and a ligand filter and performs active and inactive docks on one or more matching receptor+ligand pairs.
