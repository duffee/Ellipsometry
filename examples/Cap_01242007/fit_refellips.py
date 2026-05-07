#!/usr/bin/env python3
"""
Ellipsometry model fitting for Cap_01242007 dataset using refellips.

Layer stack (from Neha Singh's model):
    Air / Ta2O5 / EMA (Ta2O5 + void) / Ta metal

The Ta metal substrate uses point-by-point (PBP) optical constants from
ta_pbp.mat. The Ta2O5 uses a Cauchy dispersion model (transparent in
the visible/NIR). The interfacial layer between oxide and metal is an EMA
of Ta2O5 and voids, representing surface roughness from the oxidation.

Reference: Neha Singh's improved model using genosc Ta2O5 + void EMA.
The ~1 nm porous interfacial layer was found to be mostly void (>80%),
consistent with Curt's theory of columnar grain roughness from ion
mobility differences during the oxidation reaction.

Jovan Trujillo
Advanced Electronics and Photonics Core
Arizona State University
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

import refellips
from refellips import (
    DataSE,
    StructureSE,
    SlabSE,
    MixedSlabSE,
    Cauchy,
    RI,
    ReflectModelSE,
    ObjectiveSE,
)
from refnx.analysis import CurveFitter


# ============================================================
# Helper functions
# ============================================================

def load_mat_file(filepath):
    """Load a WVASE .mat file (wavelength/eV, n, k format)."""
    lines = []
    units = None
    with open(filepath, "r") as f:
        for i, line in enumerate(f):
            if i == 1:
                units = line.strip()
            if i < 3:
                continue
            parts = line.split()
            if len(parts) >= 3:
                lines.append([float(x) for x in parts[:3]])
    arr = np.array(lines)
    if units == "eV":
        wav_nm = 1239.842 / arr[:, 0]
        n = arr[:, 1]
        k = arr[:, 2]
        if wav_nm[0] > wav_nm[-1]:
            wav_nm = wav_nm[::-1]
            n = n[::-1]
            k = k[::-1]
    else:  # nm
        wav_nm = arr[:, 0]
        n = arr[:, 1]
        k = arr[:, 2]
    return wav_nm, n, k


def load_se_data(data_file):
    """Load VASE data file, return wavelength (nm), aoi, psi, delta arrays."""
    raw = []
    units = None
    with open(data_file, "r") as f:
        for i, line in enumerate(f):
            if i == 3:
                units = line.strip()
            if i < 4:
                continue
            parts = line.split()
            if not parts:
                continue
            try:
                float(parts[0])
            except ValueError:
                break
            raw.append([float(x) for x in parts[:4]])
    raw = np.array(raw)
    wavelength = raw[:, 0]
    if units == "eV":
        wavelength = 1239.842 / wavelength
    # Sort by wavelength (eV data comes in reverse order)
    sort_idx = np.argsort(wavelength)
    return wavelength[sort_idx], raw[sort_idx, 1], raw[sort_idx, 2], raw[sort_idx, 3]


def fit_wafer(wavelength, aoi, psi, delta, ta_metal_ri, wav_min, wav_max, name):
    """
    Fit a single wafer dataset with the Air/Ta2O5/EMA/Ta model.
    Returns dict with fit results and model predictions.
    """
    mask = (wavelength >= wav_min) & (wavelength <= wav_max)
    wl = wavelength[mask]
    ai = aoi[mask]
    ps = psi[mask]
    dl = delta[mask]

    data = DataSE(data=(wl, ai, ps, dl), name=name)

    # Build model
    ta2o5_cauchy = Cauchy(A=2.10, B=0.01, C=0.0, wavelength=wl, name="Ta2O5")
    air_ri = RI(dispersion=(1.0, 0.0), name="air")
    void_ri = RI(dispersion=(1.0, 0.0), name="void")

    air_layer = SlabSE(thick=0, ri=air_ri, rough=0, name="air")
    ta2o5_layer = SlabSE(thick=1960.0, ri=ta2o5_cauchy, rough=2.0, name="Ta2O5")
    ema_layer = MixedSlabSE(
        thick=1.0, ri_A=ta2o5_cauchy, ri_B=void_ri,
        vf_B=0.80, rough=0.5, name="EMA_interface",
    )
    ta_substrate = SlabSE(thick=0, ri=ta_metal_ri, rough=0.5, name="Ta_metal")

    structure = StructureSE(
        components=[air_layer, ta2o5_layer, ema_layer, ta_substrate],
        name=f"Air/Ta2O5/EMA/Ta ({name})",
        wavelength=wl,
    )

    # Fit parameters
    ta2o5_cauchy.A.setp(value=2.10, bounds=(1.8, 2.8), vary=True)
    ta2o5_cauchy.B.setp(value=0.01, bounds=(0.0, 0.1), vary=True)
    ta2o5_cauchy.C.setp(value=0.0, bounds=(0.0, 0.01), vary=True)
    ta2o5_layer.thick.setp(value=1960.0, bounds=(1500, 2500), vary=True)
    ta2o5_layer.rough.setp(value=2.0, bounds=(0, 20.0), vary=True)
    ema_layer.thick.setp(value=1.0, bounds=(0.1, 20.0), vary=True)
    ema_layer.vf_B.setp(value=0.80, bounds=(0.3, 0.99), vary=True)
    ema_layer.rough.setp(value=0.5, bounds=(0, 5.0), vary=False)
    ta_substrate.rough.setp(value=0.5, bounds=(0, 5.0), vary=False)

    model = ReflectModelSE(structure, name=name)
    objective = ObjectiveSE(model, data, use_weights=False)

    # Grid search for initial thickness
    best_chi2 = np.inf
    best_thick = 1960.0
    for thick_test in np.arange(1800, 2200, 5):
        ta2o5_layer.thick.value = thick_test
        chi2 = objective.chisqr()
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_thick = thick_test
    ta2o5_layer.thick.value = best_thick

    # LM refinement
    fitter = CurveFitter(objective)
    fitter.fit(method="least_squares")

    # Calculate MSE
    n_data = len(wl) * 2
    n_params = len(objective.varying_parameters())
    chi2 = objective.chisqr()
    reduced_chi2 = chi2 / (n_data - n_params)
    mse = np.sqrt(reduced_chi2)

    # Generate model predictions
    wav_aoi = np.column_stack([wl, ai])
    psi_model, delta_model = model(wav_aoi)

    results = {
        "name": name,
        "wavelength": wl,
        "aoi": ai,
        "psi_data": ps,
        "delta_data": dl,
        "psi_model": psi_model,
        "delta_model": delta_model,
        "chi2": chi2,
        "reduced_chi2": reduced_chi2,
        "mse": mse,
        "ta2o5_thick": ta2o5_layer.thick.value,
        "ta2o5_rough": ta2o5_layer.rough.value,
        "ema_thick": ema_layer.thick.value,
        "ema_vf": ema_layer.vf_B.value,
        "cauchy_A": ta2o5_cauchy.A.value,
        "cauchy_B": ta2o5_cauchy.B.value,
        "cauchy_C": ta2o5_cauchy.C.value,
    }
    return results


def plot_wafer_fit(results, output_file):
    """Plot Psi and Delta fit for a single wafer (one figure per angle)."""
    wl = results["wavelength"]
    ai = results["aoi"]
    angles = np.unique(ai)

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

    for idx, angle in enumerate(angles):
        mask = np.isclose(ai, angle, atol=0.1)
        wl_a = wl[mask]
        color = colors[idx % len(colors)]
        label = f"{angle:.0f}°"

        axes[0].plot(wl_a, results["psi_data"][mask], ".", color=color,
                     markersize=1.5, alpha=0.4)
        axes[0].plot(wl_a, results["psi_model"][mask], "-", color=color,
                     linewidth=1.0, label=label)

        axes[1].plot(wl_a, results["delta_data"][mask], ".", color=color,
                     markersize=1.5, alpha=0.4)
        axes[1].plot(wl_a, results["delta_model"][mask], "-", color=color,
                     linewidth=1.0, label=label)

    axes[0].set_ylabel(r"$\Psi$ (degrees)")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].set_title(f"{results['name']} — MSE = {results['mse']:.4f}")

    axes[1].set_ylabel(r"$\Delta$ (degrees)")
    axes[1].set_xlabel("Wavelength (nm)")
    axes[1].legend(loc="best", fontsize=8)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot: {output_file}")


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    # Load Ta metal optical constants (point-by-point)
    ta_metal_file = "../../data/Metal_Oxides/Ta metal/ta_pbp.mat"
    ta_wav, ta_n, ta_k = load_mat_file(ta_metal_file)
    print(f"Ta metal: {len(ta_wav)} pts, {ta_wav[0]:.1f}-{ta_wav[-1]:.1f} nm")

    ta_metal_ri = RI(dispersion=(ta_wav / 1000.0, ta_n, ta_k), name="Ta_metal")

    # Fit range: exclude deep UV where Ta2O5 absorbs
    wav_min = max(ta_wav.min(), 320.0)
    wav_max = ta_wav.max()
    print(f"Fit range: {wav_min:.1f} - {wav_max:.1f} nm\n")

    # Data files
    wafer_files = sorted(Path(".").glob("wafer_*_01242007.dat"))
    print(f"Found {len(wafer_files)} wafer data files\n")

    all_results = []

    for wf in wafer_files:
        wafer_name = wf.stem.replace("_01242007", "").replace("_", " ").title()
        print(f"Fitting {wafer_name} ({wf.name})...")

        wavelength, aoi, psi, delta = load_se_data(str(wf))
        print(f"  Data: {len(wavelength)} pts, "
              f"{wavelength.min():.0f}-{wavelength.max():.0f} nm, "
              f"angles: {np.unique(aoi).round(0).astype(int)}°")

        results = fit_wafer(wavelength, aoi, psi, delta,
                            ta_metal_ri, wav_min, wav_max, wafer_name)

        print(f"  MSE: {results['mse']:.4f}")
        print(f"  Ta2O5: {results['ta2o5_thick']:.1f} nm, "
              f"n(600nm)={results['cauchy_A'] + results['cauchy_B']/0.6**2 + results['cauchy_C']/0.6**4:.4f}")
        print(f"  EMA: {results['ema_thick']:.2f} nm, "
              f"{results['ema_vf']*100:.1f}% void")

        plot_file = f"fit_{wf.stem}.png"
        plot_wafer_fit(results, plot_file)
        all_results.append(results)
        print()

    # Summary table
    print("=" * 70)
    print(f"{'Wafer':<12} {'Ta2O5 (nm)':<12} {'EMA (nm)':<10} "
          f"{'Void %':<8} {'n@600nm':<10} {'MSE':<8}")
    print("-" * 70)
    for r in all_results:
        n600 = r["cauchy_A"] + r["cauchy_B"] / 0.6**2 + r["cauchy_C"] / 0.6**4
        print(f"{r['name']:<12} {r['ta2o5_thick']:<12.1f} {r['ema_thick']:<10.2f} "
              f"{r['ema_vf']*100:<8.1f} {n600:<10.4f} {r['mse']:<8.4f}")
    print("=" * 70)
    print("\nDone. All fits complete.")

