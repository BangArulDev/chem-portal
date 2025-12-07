class SectionHeader extends HTMLElement {

  constructor() {
    super();
  }

  connectedCallback() {
    this.render();
  }

  static get observedAttributes() {
    return ['title', 'link-text', 'link-href'];
  }

  attributeChangedCallback(name, oldValue, newValue) {
    this.render();
  }

  render() {
    this.innerHTML = '';

    const title = this.getAttribute('title') || 'Default Title';
    const linkText = this.getAttribute('link-text') || 'See All';
    const linkHref = this.getAttribute('link-href') || '#';

    const template = document.createElement('template');
    template.innerHTML = `
      <div class="section-header">
        <h3 class="section-title">${title}</h3>
        ${
          linkText
            ? `<a class="see-all-link" href="${linkHref}">${linkText}</a>`
            : ''
        }
      </div>
    `;

    this.appendChild(template.content.cloneNode(true));
  }
}

customElements.define('section-header', SectionHeader);