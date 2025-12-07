document.addEventListener("DOMContentLoaded", () => {
  let BASE_BACKEND_URL = "http://127.0.0.1:8000";

  const wait = document.getElementById("wait");
  const err_msg = document.getElementById("err-msg");
  const login_form = document.getElementById("login-form");

  fetch("/config")
    .then((resp) => resp.json())
    .then((json) => {
      const api = json.backendUrl;
      if (!api) return; // Gatau mau ngapain lagi
      BASE_BACKEND_URL = api;

      const url = new URL(api);
      url.pathname = "/profile";
      fetch(url, { headers: { Authorization: `Bearer ${localStorage.getItem("chem-token")}` } })
        .then((resp) => resp.json())
        .then((json) => {
          if (json.user !== null) {
            window.location.href = "/src/pages/mainMenu.html";
          } else {
            wait.hidden = true;
            wait.ariaHidden = "true";
            login_form.hidden = false;
            login_form.ariaHidden = "false";
          }
        })
        .catch((e) => (err_msg.innerText = e.message));
    });

  document.getElementById("show-pw").addEventListener("click", () => {
    const pw = document.getElementById("pw");
    const sp = document.getElementById("show-pw");

    if (pw.type === "password") {
      sp.src = "/src/images/eye-off.svg";
      pw.type = "text";
    } else {
      sp.src = "/src/images/eye.svg";
      pw.type = "password";
    }
  });

  document.getElementById("btn-login").addEventListener("click", (ev) => {
    ev.preventDefault();

    const name = document.getElementById("name");
    const uname = document.getElementById("uname");
    const pw = document.getElementById("pw");
    const e = document.getElementById("login_err");

    e.innerHTML = "";

    const body = JSON.stringify({ name: name.value, username: uname.value, password: pw.value });
    const method = "POST";
    const headers = { "Content-Type": "application/json" };

    const url = new URL(BASE_BACKEND_URL);
    url.pathname = "/signup";
    fetch(url, { method, body, headers })
      .then(async (resp) => ({ status: resp.status, json: await resp.json() }))
      .then(({ status, json }) => {
        if (status === 201) {
          e.innerText = "Signup berhasil, mohon ditunggu";
          window.location.href = "/src/pages/login.html";
        } else if (json.error) {
          e.innerText = json.error;
        } else {
          e.innerText = "Kesalahan server, coba lagi nanti";
        }
      })
      .catch((e) => {
        e.innerText = "Kesalahan server atau klien, coba lagi nanti";
      });
  });
});
