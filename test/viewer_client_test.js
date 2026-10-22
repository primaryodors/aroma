const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const elements = Object.fromEntries(['gobtn','pubchembtn','result','brot','rflp','smiles'].map(id => [id, {
    value: '', disabled: false, classList: {add() {}, remove() {}}
}]));
let cleared = 0;
const context = {
    document: {getElementById: id => elements[id]},
    window: {location: {pathname: '/smiles_view.php'}},
    URLSearchParams, AbortController, Blob, setTimeout, clearTimeout,
    stage: {removeAllComponents() { cleared++; }, setSpin() {},
        loadFile: async () => ({addRepresentation() {}, autoView() {}})},
    color_scheme: 'test',
    fetch: async () => ({ok: false, status: 500, text: async () => ''})
};
vm.createContext(context);
vm.runInContext(fs.readFileSync('www/js/smiles_view.js', 'utf8'), context);
(async () => {
    await context.go('CCO');
    assert.match(elements.result.value, /HTTP 500/);
    assert.equal(elements.gobtn.disabled, false);
    assert.equal(elements.pubchembtn.disabled, false);
    assert.ok(cleared >= 2);
    context.fetch = async () => ({ok: true, status: 200,
        text: async () => JSON.stringify({sdf: 'example\nM  END\n', log: 'fresh molecule'})});
    await context.go('CCO');
    assert.equal(elements.result.value, 'fresh molecule');
    context.stage.loadFile = async () => { throw new Error('Invalid structure'); };
    await context.go('CCO');
    assert.equal(elements.gobtn.disabled, false);
    assert.match(elements.result.value, /Invalid structure/);
    context.fetch = async () => { throw Object.assign(new Error('aborted'), {name: 'AbortError'}); };
    await context.go('CCO');
    assert.match(elements.result.value, /timed out/);
    assert.equal(elements.gobtn.disabled, false);
    context.stage.loadFile = async () => ({addRepresentation() {}, autoView() {}});
    context.fetch = async () => ({ok: true, text: async () => 'example\nM  END\n'});
    await context.go('CCO', 'pubchem');
    assert.match(elements.result.value, /PubChem/);
    console.log('Viewer success, HTTP failure, and display failure passed.');
})().catch(error => { console.error(error); process.exitCode = 1; });
