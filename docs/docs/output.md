# Output in MetalDock

MetalDock generates a dedicated directory named "output," which contains four distinct subdirectories, each serving a specific purpose in the docking process.

## 1. file_prep

The "file_prep" directory plays a crucial role in preparing the necessary files for your docking experiments. Within this directory, MetalDock performs various tasks, including protonation of proteins and other essential file modifications to ensure accurate docking.

## 2. QM (Quantum Mechanical Calculation)

Within the "QM" directory, MetalDock carries out Quantum Mechanical (QM) calculations, depending on your chosen settings. For each calculation, a separate subdirectory is created. This subdirectory houses either a geometry optimization or a single-point calculation. It is important to note that for cases where system convergence becomes an issue, a detailed examination of the subdirectory outputs can help in troubleshooting any problems that may arise.

## 3. docking

The "docking" directory is where the core docking process takes place. MetalDock performs the entire docking procedure here, providing all the necessary files required by the AutoDock engine. For a more detailed insight into the docking process, you can access the ".dlg" file, which contains comprehensive information about the docking run.

## 4. results

The "results" directory holds the output of the docking procedure. It presents the obtained poses in two file formats: ".pdbqt" and ".xyz." It's worth noting that AutoDock generates these structures without non-polar hydrogens. To complete the structures, you may need to manually add the missing non-polar hydrogens at a later stage.

The organization of the "output" directory ensures that you have easy access to all the essential files and results generated during the docking process, facilitating further analysis and exploration.