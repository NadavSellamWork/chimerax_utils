import os
import json
import requests
import json
from typing import List, Optional
import random
from contextlib import contextmanager
import shutil


headers = {
    "content-type": "multipart/form-data; boundary=----WebKitFormBoundary57TtBhHy2dVJKw3H",
}
# Define the form-data body
body = (
    "------WebKitFormBoundary57TtBhHy2dVJKw3H\r\n"
    "Content-Disposition: form-data; name=\"command\"\r\n\r\n"
    "COMMAND"
    "------WebKitFormBoundary57TtBhHy2dVJKw3H--\r\n"
)

class ChimeraCommandManager:
    def __init__(self, chimera_port=8889):
        print(f'Please run <chimerax --cmd "remotecontrol rest start port {chimera_port} json true"> in your command line')
        self.chimera_url = f"http://127.0.0.1:{chimera_port}/run"
        self.index = 1 
    
    def __call__(self, command):
        new_body = body.replace("COMMAND", command + "\r\n")
        output = requests.post(self.chimera_url, headers=headers, data=new_body)
        return output.json()
    
    @contextmanager
    def command_dump(self):
        """
            this context manager makes all the commands that are given within it to run as a single commmand
            where the sub commands are seperated by a semicolon
        """
        command_queue = []
        old_command_function = self.__call__
        self.__call__ = lambda command: command_queue.append(command)
        try:
            yield None
        finally:
            self.__call__ = old_command_function
            command = ";".join(command_queue)
            self(command)
    
    def get_index(self):
        return_value = self.index
        self.index += 1
        return return_value
    
    def clear(self):
        self("close all")
        self.index = 1

    def background_color(self, color="#FFFFF0"):
        self(f"set bgColor {color}")
        
    def look_at(self, objects):
        indexes = [o.index for o in objects]
        command = " ".join(["view"] + [f"#{i}" for i in indexes] )
        self(command)
    
    def load_protein(self, file_path):
        file_path = os.path.abspath(file_path)
        assert file_path.endswith(".pdb") or file_path.endswith(".cif")
        index = self.get_index()
        self(f"open {file_path}")
        protein = Protein(self, index)
        return protein
    
    def load_protein_folder(self, folder_path):
        """
            this function will load all the proteins in the folder and align them all to the first one
            returns a list of the proteins
        """
        folder_path = os.path.abspath(folder_path)
        folder_files = os.listdir(folder_path)
        protein_files = [file for file in folder_files if file.endswith(".pdb") or file.endswith(".cif")]
        protein_files = sorted(protein_files)
        proteins = [self.load_protein(os.path.join(folder_path, file)) for file in protein_files]
        for p in proteins[1:]:
            p.align(proteins[0])
        return proteins
    
    def load_and_fix_protein(self, protein_path):
        """
            this function will load a protein, fix side chains, add terminal atoms and save that protein to the same path.
            the function will return the protein object after correction
        """
        p = self.load_protein(protein_path)
        p.fix_missing_sidechains()
        p.add_terminal_atoms()
        p.save(protein_path)
        return p

    def load_density(self, file_path):
        file_path = os.path.abspath(file_path)
        assert file_path.endswith(".ccp4") or file_path.endswith(".map")
        index = self.get_index()
        self(f"open {file_path}")
        return Density(self, index)
    
    def silhouettes(self):
        self("graphics silhouettes true")
    
    def show_side_view(self):
        self.c('ui tool show "Side View"')
    
    def save_image(self,save_path:str, targets=None,transparent_background=True):
        command = "save "
        if targets is not None:
            command = command + " ".join([str(target.index) for target in targets]) +" "
        command = command + save_path + " "
        if transparent_background:
            command = command + "transparentBackground true"
        self(command)

    def create_protein_structure_series(self, proteins):
        index = self.get_index()
        proteins_morph = ProteinStructureMorph(self, index, proteins)
        return proteins_morph
    
    def create_density_series(self, densities=None,density_files=None):
        index = self.get_index()
        proteins_morph = DensityMapSeries(self, index, densities,density_files)
        return proteins_morph
    
    def load_time_series_folder(self, folder_path, align_proteins_to=None, step_size=1):
        """
            this function will accept the path of a folder which contains a time series
            this folder will have subfolders which are ordered in the frame order
            each subfolder will contains the same filenames 
            each subfolder will represent a state, made from atomic models and density maps

            align_proteins_to: a str or a protein object

            returns:
                - protein_series_objects - a list of protein_series objects
                - protein_series_names - a list of the same length, with the corepsonding files names of the protein series (from which files the protein series were created)
                - density_series_objects - a list of density series objects
                - density_series_namesa - list of the same length, with the corepsonding files names of the density series (from which files the density series were created)
                _ scene play function - a function to play and control the entire scene all together

            the scene play function will accept a start and end index, and a pause frames argument, and will play the entire scene together 
            according to these arguments
            the scene play function has the function signature of: play_scene_function(start=1,end=-1, pause_frames=10)

            an example of the folder sturcture that this function will accept: 
            - input_folder_path
                - frame_0
                    - protein_1.cif
                    - protein_2.pdb
                    - density_1.ccp4
                    - density_2.map
                - frame_1
                    - protein_1.cif
                    - protein_2.pdb
                    - density_1.ccp4
                    - density_2.map
                - frame_2
                    - protein_1.cif
                    - protein_2.pdb
                    - density_1.ccp4
                    - density_2.map

                .
                .
                .

                - frame_100
                    - protein_1.cif
                    - protein_2.pdb
                    - density_1.ccp4
                    - density_2.map
        """
        def numeric_key_function(name):
            """
                this function will be used when sorting file names, to order them based on the numbers that the basename contains
            """
            return int("".join(char for char in os.path.basename(name) if char.isnumeric()))
        subfolders = sorted(os.listdir(folder_path), key=numeric_key_function)
        subfolders = [os.path.join(folder_path, sub_folder) for sub_folder in subfolders]
        subfolders = subfolders[::step_size]
        folder_files = set(os.listdir(subfolders[0]))
        for sub_folder in subfolders[1:]:
            if not set(os.listdir(sub_folder)) == folder_files:
                raise NameError(f"folder {sub_folder} had a different content compared to {subfolders[0]}")
        trajectory_files = {file_name: [] for file_name in folder_files}
        for sub_folder in subfolders:
            for file in folder_files:
                trajectory_files[file].append(os.path.join(sub_folder, file))
        protein_trajectories = []
        protein_trajectories_names = []
        density_trajectories = []
        density_trajectory_names = []
        
        if align_proteins_to is not None:
            if isinstance(align_proteins_to, str):
                protein_to_align_to = self.load_protein(align_proteins_to)
            else:
                protein_to_align_to = align_proteins_to

        for file in folder_files:
            trajectory = trajectory_files[file]
            if file.endswith(".pdb") or file.endswith(".cif"):
                protein_trajectories_names.append(file)
                proteins = [self.load_protein(trajectory_file) for trajectory_file in trajectory]
                if align_proteins_to is not None:
                    for protein in proteins:
                        protein.align(protein_to_align_to)
                protein_series = self.create_protein_structure_series(proteins)
                protein_trajectories.append(protein_series)
            else:
                density_trajectory_names.append(file)
                density_series = self.create_density_series(density_files=trajectory)
                density_series.transparency(0.5)
                density_trajectories.append(density_series)
        if align_proteins_to is not None:
            protein_to_align_to.hide()

        def play_scene_function(start=1,end=-1, pause_frames=10):
            with self.command_dump():
                for protein_trajectory in protein_trajectories:
                    protein_trajectory.play_section(start, end, pause_frames)
                for density_trajectory in density_trajectories:
                    density_trajectory.play_section(start, end, pause_frames)
        return protein_trajectories, protein_trajectories_names, density_trajectories, density_trajectory_names, play_scene_function
        
class ChimeraObject:
    def __init__(self, c:ChimeraCommandManager, index: int):
        self.c = c 
        self.index = index
    
    def show(self):
        self.c(f"show #{self.index}")

class Protein(ChimeraObject):
    def hide(self, residue_range=None):
        # if residue_range is not None:
        #     self.c(f"hide #{self.index}:{residue_range[0]}-{residue_range[1]} atoms")
        #     self.c(f"hide #{self.index}:{residue_range[0]}-{residue_range[1]} cartoon")
        # else:
        #     self.c(f"hide #{self.index} atoms")
        #     self.c(f"hide #{self.index} cartoon")
        if residue_range is not None:
            self.c(f"hide #{self.index}:{residue_range[0]}-{residue_range[1]} models")
        else:
            self.c(f"hide #{self.index} models")
    
    def show_atoms(self, residue_range=None, backbone=None):
        command = f"show #{self.index}"
        if residue_range is not None:
           command = command + f":{residue_range[0]}-{residue_range[1]}" 
        if backbone:
            command = command + " @ca,n,c,o"
        self.c(command)
    
    def show_cartoon(self, residue_range=None):
        command = f"show #{self.index} "
        if residue_range is not None:
           command = command + f":{residue_range[0]}-{residue_range[1]}" 
        command = command + "cartoon"
        self.c(command)
    
    def color(self,color, residue_range=None):
        command = f"color #{self.index}"
        if residue_range is not None:
            command = command + f":{residue_range[0]}-{residue_range[1]}"
        command = command +f" {color}"
        self.c(command)
    
    def color_by_element(self, residue_range=None):
        command = f"color #{self.index} "
        if residue_range is not None:
            command = command + f":{residue_range[0]}-{residue_range[1]} "
        command = command + "byelement atoms"
        self.c(command)
    
    def fix_residue_indexing(self):
        self.c(f"renumber #{self.index}")
    
    def add_terminal_atoms(self):
        self.c(f"addh #{self.index}")
        self.c(f"delete #{self.index}@H*")
    
    def fix_missing_sidechains(self):
        command = f"dockprep #{self.index} ah false ac false delAltLocs false"
        self.c(command)
    
    def slice_residues(self, residue_range):
        self.c(f"split #{self.index} atoms :{residue_range[0]}-{residue_range[1]}")
        self.c(f"combine #{self.index}.1")
        return Protein(self.c, self.c.get_index())

    def align(self, other):
        self.c(f"match #{self.index} to #{other.index}")
    
    def save(self, path):
        command = f"save {path} #{self.index}"
        self.c(command)
    
    def atom_count(self):
        output = self.c(f"info #{self.index}")["json values"][0]
        d = json.loads(output)[0]
        return d["num atoms"]
    
    def color_residues(self, colors):
        output = self.c(f"info #{self.index}")["json values"][0]
        d = json.loads(output)[0]
        assert d["num residues"] == len(colors)
        
        with self.c.command_dump():
            for i, c in enumerate(colors):
                red = int(255 * (1 - c))
                blue = int(255 * c)
                green = 0
                command = f"color #{self.index}:{i+1} #{red:02X}{green:02X}{blue:02X}"
                self.c(command)
        

class Density(ChimeraObject):    
    def hide(self):
        self.c(f"volume #{self.index} hide")
    
    def level_set(self, level_set):
        self.c(f"volume #{self.index} level {level_set}")
    
    def transparency(self, transparency):
        self.c(f"volume #{self.index} transparency {transparency}")
    
    def mesh(self):
        command = f"volume #{self.index} style mesh"
        self.c(command)
    
    def surface(self):
        command = f"volume #{self.index} style surface"
        self.c(command)
    
    def color(self,color):
        command = f"volume #{self.index} color {color}"
        self.c(command)
    
    def create_slice_around_protines_resiudes(self, residue_range, proteins,backbone=False, padding=1.0, level_set=0.3, transparency=50):
        command = f"volume zone #{self.index} near "
        for protein in proteins:
            command = command + f"#{protein.index}:{residue_range[0]}-{residue_range[1]} "
            if backbone:
                command = command + "@C,CA,N,O "
        command = command + f"range {padding} newMap true"
        self.c(command)
        new_density = Density(self.c, self.c.get_index())
        new_density.level_set(level_set)
        return new_density

    def carve_density_around_protein(self, protein, padding=1.0):
        """
        Take protein structure object and carve density around it with a specific padding
        This avoids density to be around the edge
        """
        command = f"vop cover #{self.index} atom #{protein.index} pad {padding}"
        self.c(command)
        new_density = Density(self.c, self.c.get_index())
        return new_density

    def save_density(self, path):
        command = f"save {path} models #{self.index}"
        self.c(command)

class ProteinStructureMorph(Protein):
    def __init__(self, c:ChimeraCommandManager, index: int, proteins: List[Protein]):
        super().__init__(c, index)
        command = "morph " + " ".join([f"#{protein.index}" for protein in proteins]) + f" frames 1"
        self.c(command)
        self.length = len(proteins)
    
    def set_frame(self, frame_index: int):
        command = f"coordset #{self.index} {frame_index}"
        self.c(command)
    
    def play_section(self, start=1,end=-1, pause_frames=10):
        command = f"coordset #{self.index} {start}"
        if end:
            command += f",{end}"
        if pause_frames is not None and pause_frames > 0:
            command += f" pauseFrames {pause_frames}"
        self.c(command)
    
    def __len__(self):
        return self.length
    
    def __getitem__(self, sub_index):
        assert isinstance(sub_index, int) and 1 <= sub_index <= self.length
        index = f"{self.index}.{sub_index}"
        protein = Protein(self.c, index)
        return protein

class DensityMapSeries(Density):
    def __init__(self, c:ChimeraCommandManager, index: int, densities: Optional[List[Density]],density_files:Optional[List[str]]):
        # make sure that only one of them is not None
        assert ((densities and density_files) is None) and ((densities or density_files) is not None)
        super().__init__(c, index)
        random_number = random.random()
        if densities is not None:
            self.save_path = os.path.abspath(f"/tmp/temp_{random_number:.6f}")
            os.makedirs(self.save_path)
            # for a very stupid reason, chimerax cannot create a vseries from existing volumes, and they must be loaded using the open command
            save_command = f"save {os.path.join(self.save_path, '%d.ccp4')} models #{','.join([str(density.index) for density in densities])}"
            self.c(save_command)
            file_paths = [os.path.join(self.save_path, file_name) for file_name in os.listdir(self.save_path)]
            file_paths = sorted(file_paths, key=lambda file_path: int(os.path.basename(file_path).split(".")[0]))
            open_command = f'open {" ".join(file_paths)} vseries true'
            self.c(open_command)
            self.length = len(densities)
        else:
            open_command = f'open {" ".join(density_files)} vseries true'
            self.c(open_command)
            self.length = len(density_files)
            
        self.show_all_timestep()
        if densities is not None:
            # after the show_all_timestep function is done, all the models were loaded into the cache, and we can delete the temporary files
            shutil.rmtree(self.save_path)
    
    def show_all_timestep(self):
        indexes = [f"#{self.index}.{i}" for i in range(1, self.length + 1)]
        command = f"volume {' '.join(indexes)} show"
        self.c(command)
        command = f"volume {' '.join(indexes)} hide"
        self.c(command)
        self.set_frame(0)
    
    def set_frame(self, frame_index: int):
       command = f"vseries play #{self.index} jumpTo {frame_index}"
       self.c(command)
    
    def play_section(self, start=1,end=-1, pause_frames=10):
        if end < 0:
            end = end + 1 + self.length
        command = f"vseries play #{self.index} range {start},{end}"
        if pause_frames is not None and pause_frames > 0:
            command = command + f" pauseFrames {pause_frames}"
        return self.c(command)
    
    def __len__(self):
        return self.length
    
    def __getitem__(self, sub_index):
        assert isinstance(sub_index, int) and 1 <= sub_index <= self.length
        index = f"{self.index}.{sub_index}"
        density = Density(self.c, index)
        return density
