document.querySelectorAll('div.sourceCode').forEach(block => {
  const code = block.querySelector('code');
  if (!code) return;
  const button = document.createElement('button');
  button.className = 'copy-code'; button.type = 'button'; button.textContent = 'Copy';
  button.setAttribute('aria-label', 'Copy code to clipboard');
  button.addEventListener('click', async () => {
    try { await navigator.clipboard.writeText(code.innerText); button.textContent = 'Copied'; }
    catch { button.textContent = 'Select code to copy'; }
    setTimeout(() => { button.textContent = 'Copy'; }, 2000);
  });
  block.append(button);
});
const search = document.querySelector('#course-search');
if (search) {
  const results = document.querySelector('#search-results');
  const status = document.querySelector('#search-status');
  let pages;
  async function update() {
    const query = search.value.trim().toLowerCase();
    results.replaceChildren();
    if (!query) { status.textContent = 'Search lecture notes, assignments, and course resources.'; return; }
    try {
      if (!pages) { const response = await fetch('search-index.json'); if (!response.ok) throw Error(); pages = await response.json(); }
      if (query !== search.value.trim().toLowerCase()) return;
      const terms = query.split(/\s+/);
      const matches = pages.filter(page => terms.every(term => (page.title + ' ' + page.text).toLowerCase().includes(term)));
      status.textContent = matches.length + (matches.length === 1 ? ' result' : ' results');
      for (const page of matches) {
        const article = document.createElement('article'); article.className = 'search-result';
        const heading = document.createElement('h2'); const link = document.createElement('a');
        link.href = page.url; link.textContent = page.title; heading.append(link);
        const snippet = document.createElement('p');
        const start = Math.max(0, page.text.toLowerCase().indexOf(terms[0]) - 60);
        snippet.textContent = (start ? '…' : '') + page.text.slice(start, start + 230) + '…';
        article.append(heading, snippet); results.append(article);
      }
    } catch { status.textContent = 'Search could not load. Please use the Lectures or Assignments links.'; }
  }
  search.addEventListener('input', update);
}
