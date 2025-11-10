http://www.primaryodors.org


# DEPENDENCIES

g++ openbabel python3-pybel librdkit1 python3-rdkit rdkit-data php php-curl php-gd


# Performing Docks

First navigate to the AROMA root folder and run `make`. Note all command line tools must be run from this project root folder.

Next, type `./dock.sh` with the receptor as the first argument and the odorant as the second argument.

In lieu of a receptor, you may specify:
- `*` meaning to run for all receptors;
- `emp` meaning all empirical pairs for the odorant(s);
- `ago` meaning all receptors where each odorant is a known agonist;
- `e1` ~ `e100` meaning only receptors whose rate of expression in the human population is known to be at least the number (percentage) specified.

Regarding expression rates, it was found by Verbeurgt et al (2014) that each olfactory receptor is expressed in a subset of individuals.
They searched whole olfactory mucosa of 26 cadavers for RNA signatures of the ORs and found that nobody expresses the full set of 400+ receptors.

In lieu of an odorant, you may specify:
- `*` meaning all odorants in the database;
- `emp` meaning all empirical pairs for the receptor(s);
- `ago` meaning all known agonists of the receptor(s);
- a perceptual note from the database, such as `citrus`, `malty`, or `acrid`.


# Web Application

You may optionally host your own AROMA web interface.

To enable the web app:
- Either set up a local web server or checkout this repository in a folder on a web host.
- Make sure your server has the `php`, `php-curl`, `php-gd`, and `openbabel` packages installed.
- After installing `php-curl`, it's important to restart the web service e.g. `sudo apache2ctl -k restart`.
- Then open the `www/symlink.sh` file in a text editor, make sure the destination folder is correct (by default it will show `/var/www/html/`
  which is usually correct for Apache2 installations), make sure you have write permissions in the 
  folder (or use `sudo`), and execute `www/symlink.sh` in a command line.
- The `data` and `www/assets` folders and all contents must also be recursively made writable by the web user.
- If on a local server, you will now have an instance of the web app at http://127.0.0.1/primarydock/ whereas if you are using a web host
  then you may have to configure your hosting to point one of your registered domains or subfolders to the `aroma/www` folder.

If you get a 403 Forbidden error, please make sure that every containing folder of the `aroma/www` folder has public execute access.
