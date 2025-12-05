// D:\chem_portal_ui\public\script.js

document.addEventListener("DOMContentLoaded", () => {
  // --- 1. Konfigurasi dan Elemen DOM ---
  const btnFind = document.getElementById("btn-find");

  // Elemen Input
  const inputSolubility = document.getElementById("solubility");
  const inputViscosity = document.getElementById("viscosity");
  const inputThermal = document.getElementById("thermal");
  const inputBoiling = document.getElementById("boiling");

  // Elemen Output
  const resName = document.getElementById("res-name");
  const resSMILES = document.getElementById("res-smiles");
  const resStatus = document.getElementById("res-status");
  const molImg = document.getElementById("mol-img");
  const downloadMolLink = document.getElementById("download-mol");
  const resProps = document.getElementById("res-props");
  const resJust = document.getElementById("res-just");

  let API_URL = "http://127.0.0.1:8000/discover";
  let BASE_BACKEND_URL = "http://127.0.0.1:8000";
  const DEFAULT_IMG_SRC =
    "data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='220' height='160'><rect width='100%' height='100%' fill='%23f6f9fb'/><text x='50%' y='50%' dominant-baseline='middle' text-anchor='middle' fill='%238a8f96' font-family='sans-serif' font-size='14'>Molecule preview</text></svg>";

  // --- 0. Load Config ---
  async function loadConfig() {
    try {
      const response = await fetch("/config");
      const config = await response.json();
      if (config.backendUrl) {
        BASE_BACKEND_URL = config.backendUrl;
        API_URL = `${BASE_BACKEND_URL}/discover`;
        console.log("Config loaded:", config);
      }
    } catch (error) {
      console.warn("Gagal memuat config, menggunakan default localhost.");
    }
  }

  // Panggil loadConfig segera
  loadConfig();

  // --- 2. Fungsi Reset UI ---
  function resetResultsUI(statusText = "Siap menerima kriteria...") {
    resName.textContent = "-";
    resSMILES.textContent = "-";
    resStatus.textContent = statusText;
    molImg.src = DEFAULT_IMG_SRC;
    downloadMolLink.href = "#";
    downloadMolLink.setAttribute("aria-disabled", "true");
    resProps.textContent = "Properti akan terdaftar di sini setelah pencarian.";
    resJust.textContent = "Hasil justifikasi akan muncul di sini.";
  }

  // --- 3. Fungsi Tampil Hasil ---
  function displayResults(data) {
    const compound = data.recommended_compound;
    const props = compound.Properties;

    // 1. Informasi Utama
    resName.textContent = props.IUPAC || "Senyawa Baru";
    resSMILES.textContent = compound.SMILES;
    resStatus.textContent = "Rekomendasi berhasil!";

    // 2. Visualisasi
    const imagePath2D = compound.Structure_2D_Path
      ? `${BASE_BACKEND_URL}${compound.Structure_2D_Path}`
      : DEFAULT_IMG_SRC;
    const imagePath3D = compound.Structure_3D_Path
      ? `${BASE_BACKEND_URL}${compound.Structure_3D_Path}`
      : "#";

    molImg.src = imagePath2D;
    molImg.alt = `Struktur 2D ${compound.SMILES}`;

    if (imagePath3D !== "#") {
      downloadMolLink.href = imagePath3D;
      downloadMolLink.setAttribute("aria-disabled", "false");
    } else {
      downloadMolLink.setAttribute("aria-disabled", "true");
    }

    // 3. Properti (Dibuat menjadi string format JSON/Preformatted)
    const propsToDisplay = {
      IUPAC: props.IUPAC,
      SMILES: props.SMILES,
      "Mol. Weight": props.MolWeight.toFixed(3),
      "Exact Mass": props.ExactMass.toFixed(3),
      "H-Bond Donors": props.HBondDonors,
      "H-Bond Acceptors": props.HBondAcceptors,
      TPSA: props.TPSA.toFixed(2),
      "Log P": props.LogP.toFixed(3),
      "Rotatable Bonds": props.RotatableBonds,
      InChI: props.InChI,
      InChIKey: props.InChIKey,
    };

    resProps.textContent = JSON.stringify(propsToDisplay, null, 2);

    // 4. Justifikasi AI
    resJust.textContent = data.justification_ai;
  }

  // --- 4. Event Listener Utama ---
  btnFind.addEventListener("click", async () => {
    resetResultsUI("Mencari senyawa baru...");

    // Kumpulkan data dari input
    const data = {
      Solubility: parseFloat(inputSolubility.value),
      Viscosity_cP: parseFloat(inputViscosity.value),
      ThermalStability_Score: parseFloat(inputThermal.value),
      BoilingPoint_C: parseFloat(inputBoiling.value),
    };

    if (Object.values(data).some(isNaN)) {
      resetResultsUI("ERROR: Semua input harus diisi dengan angka valid.");
      return;
    }

    try {
      const response = await fetch(API_URL, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(data),
      });

      const result = await response.json();

      if (!response.ok) {
        // Tangani error dari FastAPI (misalnya HTTPException)
        const errorDetail =
          result.detail || "Terjadi kesalahan tidak diketahui di backend.";
        resetResultsUI(`ERROR: ${response.status} - ${errorDetail}`);
        return;
      }

      displayResults(result);
    } catch (error) {
      console.error("Fetch error:", error);
      resetResultsUI(
        `FATAL ERROR: Gagal terhubung ke API (Pastikan FastAPI berjalan di ${BASE_BACKEND_URL}).`
      );
    }
  });

  // Panggil reset saat startup
  resetResultsUI();
});
