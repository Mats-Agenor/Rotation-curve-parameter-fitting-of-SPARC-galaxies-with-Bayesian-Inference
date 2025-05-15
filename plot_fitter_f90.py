import os
import corner
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# =============================================
# Configurações de estilo para os gráficos
# =============================================
plt.style.use('default')
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams['font.size'] = 12
plt.rc('font', **{'family':'serif', 'serif':['Times']})
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['text.usetex'] = False
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['xtick.major.width'] = 0.79
mpl.rcParams['xtick.minor.width'] = 0.79
mpl.rcParams['ytick.major.width'] = 0.79
mpl.rcParams['ytick.minor.width'] = 0.79
#plt.rcParams['figure.constrained_layout.use'] = True


def read_fortran_results(file_path):
    results = {
        "params": None,
        "rotation_curve": None,
        "profile": None
    }

    with open(file_path, "r") as file:
        lines = file.readlines()

    params = {}
    rotation_data = []
    profile = None

    for line in lines:
        line = line.strip()
        
        # Skip empty lines
        if not line:
            continue
            
        # Get profile type
        if line.startswith("# Perfil:"):
            results["profile"] = line.split(":")[1].strip()
            continue
            
        # Read parameters
        if "rs =" in line:
            parts = [p.strip() for p in line.split("±")]
            params["rs"] = float(parts[0].split("=")[1].strip())
            params["rs_std"] = float(parts[1].strip())
            continue
            
        if "rho_s =" in line:
            parts = [p.strip() for p in line.split("±")]
            val_part = parts[0].split("=")[1].strip().replace("E", "e")
            std_part = parts[1].strip().replace("E", "e")
            params["rho_s"] = float(val_part)
            params["rho_s_std"] = float(std_part)
            continue
            
        if "yd =" in line:
            parts = [p.strip() for p in line.split("±")]
            params["yd"] = float(parts[0].split("=")[1].strip())
            params["yd_std"] = float(parts[1].strip())
            continue
            
        # Read rotation curve data
        if line.startswith("# R[kpc]"):
            # The next lines should be data
            data_started = False
            continue
            
        if any(c.isdigit() for c in line) and not line.startswith("#"):
            try:
                values = list(map(float, line.split()))
                if len(values) == 7:
                    rotation_data.append(values)
            except ValueError:
                continue

    if params:
        results['params'] = params

    if rotation_data:
        rotation_data = np.array(rotation_data)
        results['rotation_curve'] = (
            rotation_data[:, 0],  # R
            rotation_data[:, 1],  # Vobs
            rotation_data[:, 2],  # e_Vobs
            rotation_data[:, 3],  # Vgas
            rotation_data[:, 4],  # Vdisk
            rotation_data[:, 5],  # Vhalo
            rotation_data[:, 6],  # Vtotal
        )

    return results


def plot_rotation_curve(file_name, results, output_dir):
    if not results or 'rotation_curve' not in results or results['rotation_curve'] is None:
        print(f"  [Aviso] Dados de curva de rotação ausentes para {file_name}")
        return

    R, Vobs, e_Vobs, Vgas, Vdisk, Vhalo, Vtotal = results['rotation_curve']
    profile = results.get('profile', file_name.split('_')[0])

    plt.figure(figsize=(8, 6))
    plt.errorbar(R, Vobs, yerr=e_Vobs, fmt='o', label="Observado", color='black')
    plt.plot(R, Vgas, label="Gás", linestyle='--')
    plt.plot(R, Vdisk, label="Disco", linestyle='--')
    plt.plot(R, Vhalo, label="Halo", linestyle='--')
    plt.plot(R, Vtotal, label="Total", linestyle='-')

    plt.xlabel("Raio (kpc)")
    plt.ylabel("Velocidade (km/s)")
    plt.title(f"Curva de Rotação - {profile}")
    plt.legend()

    plt.tight_layout()

    output_path = os.path.join(output_dir, f"{file_name}_rotation_curve.png")
    plt.savefig(output_path)
    plt.close()


def generate_pdfs_from_params(results):
    if not results or 'params' not in results or results['params'] is None:
        return None

    params = results['params']
    
    # Get parameters with defaults
    yd = params.get('yd', 0)
    yd_std = max(params.get('yd_std', 1e-6), 1e-6)
    rs = params.get('rs', 0)
    rs_std = max(params.get('rs_std', 1e-6), 1e-6)
    rho_s = params.get('rho_s', 1e-6)
    rho_s_std = max(params.get('rho_s_std', 1e-6), 1e-6)

    # Create ranges for PDFs
    yd_values = np.linspace(yd - 3*yd_std, yd + 3*yd_std, 100)
    rs_values = np.linspace(rs - 3*rs_std, rs + 3*rs_std, 100)
    rho_s_values = np.linspace(rho_s - 3*rho_s_std, rho_s + 3*rho_s_std, 100)

    # Calculate PDFs
    pdf_yd = np.exp(-0.5 * ((yd_values - yd) / yd_std) ** 2)
    pdf_rs = np.exp(-0.5 * ((rs_values - rs) / rs_std) ** 2)
    pdf_rho_s = np.exp(-0.5 * ((rho_s_values - rho_s) / rho_s_std) ** 2)

    # Normalize PDFs
    pdf_yd /= np.trapz(pdf_yd, yd_values)
    pdf_rs /= np.trapz(pdf_rs, rs_values)
    pdf_rho_s /= np.trapz(pdf_rho_s, rho_s_values)

    return {
        "yd": (yd_values, pdf_yd),
        "rs": (rs_values, pdf_rs),
        "rho_s": (rho_s_values, pdf_rho_s)
    }


def plot_pdfs(file_name, results, output_dir):
    pdfs = generate_pdfs_from_params(results)
    if not pdfs:
        print(f"  [Aviso] PDFs ausentes para {file_name}")
        return

    profile = results.get('profile', file_name.split('_')[0])
    params = results['params']
    
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot YD PDF with 1-sigma region
    yd, yd_std = params['yd'], params['yd_std']
    axs[0].plot(pdfs["yd"][0], pdfs["yd"][1], lw=2)
    axs[0].axvline(yd, color='k', linestyle='--', lw=1)
    axs[0].axvspan(yd - yd_std, yd + yd_std, color='gray', alpha=0.3)
    axs[0].set_title("PDF de YD")
    axs[0].set_xlabel("YD")
    axs[0].set_ylabel("Probabilidade")
    
    # Plot rs PDF with 1-sigma region
    rs, rs_std = params['rs'], params['rs_std']
    axs[1].plot(pdfs["rs"][0], pdfs["rs"][1], lw=2)
    axs[1].axvline(rs, color='k', linestyle='--', lw=1)
    axs[1].axvspan(rs - rs_std, rs + rs_std, color='gray', alpha=0.3)
    axs[1].set_title("PDF de rs")
    axs[1].set_xlabel("rs (kpc)")
    axs[1].set_ylabel("Probabilidade")
    
    # Plot rho_s PDF with 1-sigma region
    rho_s, rho_s_std = params['rho_s'], params['rho_s_std']
    axs[2].plot(pdfs["rho_s"][0], pdfs["rho_s"][1], lw=2)
    axs[2].axvline(rho_s, color='k', linestyle='--', lw=1)
    axs[2].axvspan(rho_s - rho_s_std, rho_s + rho_s_std, color='gray', alpha=0.3)
    axs[2].set_title("PDF de ρ_s")
    axs[2].set_xlabel("ρ_s (M☉/kpc³)")
    axs[2].set_ylabel("Probabilidade")

    plt.suptitle(f"Distribuições de Probabilidade - {profile}", y=1.05)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{file_name}_pdfs.png")
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

def plot_corner(file_name, results, output_dir):
    if not results or 'params' not in results or results['params'] is None:
        print(f"  [Aviso] Parâmetros ausentes para {file_name}")
        return

    params = results['params']
    profile = results.get('profile', file_name.split('_')[0])
    
    # Get parameters with defaults
    yd = params.get('yd', 0)
    yd_std = max(params.get('yd_std', 1e-6), 1e-6)
    rs = params.get('rs', 0)
    rs_std = max(params.get('rs_std', 1e-6), 1e-6)
    rho_s = max(params.get('rho_s', 1e-6), 1e-6)
    rho_s_std = max(params.get('rho_s_std', 1e-6), 1e-6)

    # Handle log scale for rho_s
    log_rho_s = np.log10(rho_s) if rho_s > 0 else -10
    log_rho_s_std = rho_s_std / (rho_s * np.log(10)) if rho_s > 0 else 1e-6

    # Generate samples
    samples = np.random.multivariate_normal(
        mean=[yd, rs, log_rho_s],
        cov=np.diag([yd_std**2, rs_std**2, log_rho_s_std**2]),
        size=5000
    )

    # Create corner plot
    fig = corner.corner(
        samples,
        labels=[r"$Y_D$", r"$r_s$ [kpc]", r"$\log_{10}(\rho_s)$ [M☉/kpc³]"],
        truths=[yd, rs, log_rho_s],
        quantiles=[0.16, 0.5, 0.84],  # Fixed: now has exactly 3 values
        show_titles=True,
        title_fmt=".2f",
        title_kwargs={"fontsize": 12},
        label_kwargs={"fontsize": 14}
    )

    fig.suptitle(f"Distribuição Conjunta - {profile}", y=1.02, fontsize=16)
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, f"{file_name}_corner.png")
    fig.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)


def process_all_results(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for file in os.listdir(input_dir):
        if file.endswith(".txt"):
            file_path = os.path.join(input_dir, file)
            file_name = os.path.splitext(file)[0]

            print(f"Processando {file}...")
            results = read_fortran_results(file_path)
            
            if not results:
                print(f"  [Erro] Falha ao processar {file}")
                continue

            plot_rotation_curve(file_name, results, output_dir)
            plot_pdfs(file_name, results, output_dir)
            plot_corner(file_name, results, output_dir)


if __name__ == "__main__":
    input_directory = "/home/matheus-agenor/Documentos/program/fitter_fortran/results"
    output_directory = "/home/matheus-agenor/Documentos/program/fitter_fortran/plots"
    process_all_results(input_directory, output_directory)
