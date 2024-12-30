from chimerax_utils import * 
import os 

movie_folder_path = "/Users/nadavsellam/Documents/research/proteinx_guidance/optimize_msa_and_latent_experiment_results"
c = ChimeraCommandManager()

def load_protein_structures_0():
    sub_folders_paths = [os.path.join(movie_folder_path, subfolder) for subfolder in sorted(os.listdir(movie_folder_path), key=lambda name: int(''.join([char for char in os.path.basename(name) if char.isdigit()])))]
    proteins = []
    for sub_folder in sub_folders_paths:
        file_name = os.path.join(sub_folder, "0.cif")
        protein = c.load_protein(file_name)
        proteins.append(protein)
    return proteins

def load_densities_0():
    sub_folders_paths = [os.path.join(movie_folder_path, subfolder) for subfolder in sorted(os.listdir(movie_folder_path), key=lambda name: int(''.join([char for char in os.path.basename(name) if char.isdigit()])))]
    densities = []
    for sub_folder in sub_folders_paths:
        file_name = os.path.join(sub_folder, "0.ccp4")
        density = c.load_density(file_name)
        density.color("blue")
        density.transparency(0.5)
        densities.append(density)
    return densities

def main():
    c.clear()
    proteins = load_protein_structures_0()
    altloc_a = c.load_protein("/Users/nadavsellam/Documents/research/proteinx_guidance/pdbs/1lu4_A_aligned.pdb")
    altloc_a.color("red")
    c.look_at([altloc_a])
    for protein in proteins:
        protein.align(altloc_a)
        protein.color("white")
    proteins_morph = c.create_protein_structure_series(proteins)
    proteins_morph.set_frame(0)
    proteins_morph.hide()
    proteins_morph.show_atoms()

    densities = load_densities_0()
    density_morph = c.create_density_series(densities, 2)
    density_morph.transparency(0.5)
    density_morph.level_set(0.6)
    for density in densities:
        density.hide()
    c.look_at([altloc_a])
    altloc_a.hide()

    with c.command_dump():
        proteins_morph.play_section(pause_frames=30)
        density_morph.play_section(pause_frames=30)

def main_2():
    c.clear()
    protein_morphs,protein_names, density_morphs, density_names, play_scene_function = c.load_time_series_folder(movie_folder_path, align_proteins_to="/Users/nadavsellam/Documents/research/proteinx_guidance/pdbs/1lu4_A_aligned.pdb")
    play_scene_function()
    a = 2


if __name__ == "__main__":
    main_2()