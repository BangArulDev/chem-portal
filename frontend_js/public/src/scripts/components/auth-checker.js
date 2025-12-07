const wait = document.getElementById("wait");
const logout = document.getElementById("logout");
const uname = document.getElementById("user_name");

logout.addEventListener("click", () => {
  localStorage.removeItem("chem-token");
  window.location.href = "/src/pages/login.html";
});

fetch("/config")
  .then((resp) => resp.json())
  .then((json) => {
    const api = json.backendUrl;
    if (!api) return; // Gatau mau ngapain lagi

    const url = new URL(api);
    url.pathname = "/profile";
    fetch(url, { headers: { Authorization: `Bearer ${localStorage.getItem("chem-token")}` } })
      .then((resp) => resp.json())
      .then((json) => {
        if (json.user === null) {
          window.location.href = "/src/pages/login.html";
        } else {
          console.log(`Logged in as ${json.user.name}`);

          wait.hidden = true;
          wait.ariaHidden = "true";
          logout.hidden = false;
          logout.ariaHidden = "false";
          uname.hidden = false;
          uname.ariaHidden = "false";
          uname.innerText = json.user.name;
        }
      })
      .catch((err) => {
        console.error(err);
      });
  });
