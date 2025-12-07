document.addEventListener("DOMContentLoaded", () => {
  const form = document.querySelector("form");
  const resultsContainer = document.getElementById("results");
  const headerStatus = document.getElementById("header-status");

  // Output Elements
  const resName = document.getElementById("res-name");
  const resSmiles = document.getElementById("res-smiles");
  const resStatus = document.getElementById("res-status");
  const molImg = document.getElementById("mol-img");
  const downloadMol = document.getElementById("download-mol");
  const resPropsContainer = document.getElementById("res-props");
  const resJustification = document.getElementById("res-just");
  const submitBtn = document.getElementById("submit-btn");

  const API_BASE_URL = ""; // This is now handled by relative paths

  // Mapping for solubility
  const solubilityMap = {
    rendah: 0.3,
    sedang: 0.6,
    tinggi: 0.9,
  };

  const setUIState = (state, message = "") => {
    switch (state) {
      case "loading":
        headerStatus.textContent = "Mencari...";
        headerStatus.className = "badge-status loading";
        submitBtn.disabled = true;
        resStatus.textContent = "Memproses permintaan...";
        break;
      case "success":
        headerStatus.textContent = "Selesai";
        headerStatus.className = "badge-status success";
        submitBtn.disabled = false;
        resStatus.textContent = "Prediksi Berhasil";
        break;
      case "error":
        headerStatus.textContent = "Gagal";
        headerStatus.className = "badge-status error";
        submitBtn.disabled = false;
        resStatus.textContent = "Terjadi Kesalahan";
        resJustification.innerHTML = `<p class="text-red-500">${message}</p>`;
        break;
      case "idle":
      default:
        headerStatus.textContent = "Menunggu Input";
        headerStatus.className = "badge-status";
        submitBtn.disabled = false;
        break;
    }
  };

  const renderProperties = (properties) => {
    resPropsContainer.innerHTML = ""; // Clear previous results
    const propsGrid = [
        { label: "Molecular Weight", value: properties.MolWeight?.toFixed(2) },
        { label: "Exact Mass", value: properties.ExactMass?.toFixed(2) },
        { label: "LogP", value: properties.LogP?.toFixed(2) },
        { label: "TPSA", value: properties.TPSA?.toFixed(2) },
        { label: "H-Bond Donors", value: properties.HBondDonors },
        { label: "H-Bond Acceptors", value: properties.HBondAcceptors },
        { label: "Rotatable Bonds", value: properties.RotatableBonds },
    ];

    propsGrid.forEach(prop => {
        const propEl = document.createElement("div");
        propEl.className = "prop-item";
        propEl.innerHTML = `
            <span class="prop-label">${prop.label}</span>
            <span class="prop-value">${prop.value ?? '-'}</span>
        `;
        resPropsContainer.appendChild(propEl);
    });
  };

  submitBtn.addEventListener("click", async (e) => {
    if (!form.checkValidity()) {
      form.reportValidity();
      return;
    }
    
    setUIState("loading");

    const boilingPoint = parseFloat(document.getElementById("boiling_point").value);
    const viscosity = parseFloat(document.getElementById("viscosity").value);
    const solubilityKey = document.getElementById("solubility").value;
    const thermalStability = parseFloat(document.getElementById("thermal_stability").value);

    if (!solubilityKey) {
        setUIState("error", "Silakan pilih tingkat kelarutan.");
        return;
    }

    const payload = {
      BoilingPoint_C: boilingPoint,
      Viscosity_cP: viscosity,
      Solubility: solubilityMap[solubilityKey],
      ThermalStability_Score: thermalStability,
    };

    try {
      const response = await fetch(`/api/discover`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(payload),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || `HTTP error! status: ${response.status}`);
      }

      const data = await response.json();
      
      // Update UI with new data
      resName.textContent = data.recommended_compound.Properties.IUPAC || "N/A";
      resSmiles.textContent = data.recommended_compound.SMILES;
      
      if (data.recommended_compound.Structure_2D_Path) {
        molImg.src = data.recommended_compound.Structure_2D_Path;
      } else {
        // Fallback image if not provided
        molImg.src = "data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='220' height='160'><rect width='100%' height='100%' fill='%23f8fafc'/><text x='50%' y='50%' dominant-baseline='middle' text-anchor='middle' fill='%2394a3b8' font-family='sans-serif' font-size='14'>Preview Not Available</text></svg>";
      }
      
      if (data.recommended_compound.Structure_3D_Path) {
        downloadMol.href = data.recommended_compound.Structure_3D_Path;
        downloadMol.classList.remove("disabled");
      } else {
        downloadMol.href = "#";
        downloadMol.classList.add("disabled");
      }

      renderProperties(data.recommended_compound.Properties);
      resJustification.innerHTML = `<p>${data.justification_ai.replace(/\n/g, '<br>')}</p>`;

      setUIState("success");
      console.log("UI Update complete. Data:", data);

    } catch (error) {
      console.error("Fetch error:", error);
      setUIState("error", `Gagal menghubungi API. Pastikan server backend berjalan di ${API_BASE_URL}. Error: ${error.message}`);
    }
  });
});
