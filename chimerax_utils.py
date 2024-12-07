import os
import json
import requests
import json
import time

PORT = 8889
print(f'Please run <chimerax --cmd "remotecontrol rest start port {PORT} json true"> in your command line')
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
    def __init__(self, chimera_port=PORT):
        self.chimera_url = f"http://127.0.0.1:{chimera_port}/run"
        self.index = 1 
    
    def __call__(self, command):
        new_body = body.replace("COMMAND", command + "\r\n")
        requests.post(self.chimera_url, headers=headers, data=new_body)
    
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
        folder_files = os.listdir(folder_path)
        protein_files = [file for file in folder_files if file.endswith(".pdb") or file.endswith(".cif")]
        protein_files = sorted(protein_files)
        proteins = [self.load_protein(os.path.join(folder_path, file)) for file in protein_files]
        for p in proteins[1:]:
            p.align(proteins[0])
        return proteins
    
    def load_density(self, file_path):
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
        
class ChimeraObject:
    def __init__(self, c:ChimeraCommandManager, index: int):
        self.c = c 
        self.index = index
    
    def show(self):
        self.c(f"show #{self.index}")

class Protein(ChimeraObject):
    def hide(self, residue_range=None):
        if residue_range is not None:
            self.c(f"hide #{self.index}:{residue_range[0]}-{residue_range[1]} atoms")
            self.c(f"hide #{self.index}:{residue_range[0]}-{residue_range[1]} cartoon")
        else:
            self.c(f"hide #{self.index} atoms")
            self.c(f"hide #{self.index} cartoon")
    
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
        command = f"dockprep #{self.index} ah false ac false "
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
    
