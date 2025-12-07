class PackageItem extends HTMLElement {
  _package = {
    id: null,
    name: null,
    details: null,
    price: 0.00,
    icon: null,
  };

  constructor() {
    super();
  }

  _emptyContent() {
    this.innerHTML = '';
  }

  set package(value) {
    this._package = value;
    this.render();
  }

  get package() {
    return this._package;
  }
  render() {
    this._emptyContent();
    this.innerHTML += `
      <div class="package-item">
        <div class="package-item-header">
          <div>
            <h4 class="package-item-title">${this._package.name}</h4>
            <p class="package-item-details">${this._package.details}</p>
          </div>
          <div class="package-item-icon-wrapper">
            <span class="material-symbols-outlined package-item-icon">${this._package.icon}</span>
          </div>
        </div>
        <div class="package-item-footer">
          <p class="package-item-price">Rp${this._package.price.toFixed(3)}</p>
          <button class="button button-small">Buy Now</button>
        </div>
      </div>
    `;
  }
}

customElements.define('package-item', PackageItem);