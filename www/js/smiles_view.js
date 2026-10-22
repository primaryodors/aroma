async function viewer_request(smiles, source) {
    const controller = new AbortController();
    const deadline = setTimeout(() => controller.abort(), 30000);
    try {
        const brot = document.getElementById('brot').value;
        const rflp = document.getElementById('rflp').value;
        if (source === 'pubchem') {
            if (smiles.includes('{') || brot.trim() || rflp.trim()) {
                throw new Error('PubChem accepts standard SMILES without local transformations.');
            }
            const response = await fetch('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/SDF?record_type=3d', {
                method: 'POST', body: new URLSearchParams({smiles}), signal: controller.signal
            });
            const sdf = await response.text();
            if (!response.ok || !sdf.includes('M  END')) {
                throw new Error('PubChem could not supply a 3D conformer for this SMILES.');
            }
            return {sdf, log: 'Loaded a 3D conformer from PubChem.'};
        }
        const response = await fetch(window.location.pathname, {
            method: 'POST', body: new URLSearchParams({smiles, brot, rflp}),
            cache: 'no-store', signal: controller.signal
        });
        const text = await response.text();
        let result;
        try { result = JSON.parse(text); }
        catch (_) { throw new Error('Conversion failed (HTTP ' + response.status + '). Check the host PHP settings or use PubChem lookup.'); }
        if (!response.ok || result.error) throw new Error(result.error || 'Molecule generation failed.');
        if (!result.sdf || !result.sdf.includes('M  END')) throw new Error('No valid molecular structure was returned.');
        return result;
    } catch (error) {
        if (error.name === 'AbortError') throw new Error('Molecule generation timed out. Try again or use PubChem lookup.');
        throw error;
    } finally {
        clearTimeout(deadline);
    }
}

async function go(smiles, source = 'local') {
    const buttons = ['gobtn', 'pubchembtn'].map(id => document.getElementById(id));
    const output = document.getElementById('result');
    buttons.forEach(button => { button.disabled = true; button.classList.add('grayed'); });
    output.value = source === 'pubchem' ? 'Looking up a 3D conformer...' : 'Generating molecule...';
    try {
        stage.removeAllComponents();
        stage.setSpin(false);
        const result = await viewer_request(smiles, source);
        const component = await stage.loadFile(new Blob([result.sdf], {type: 'text/plain'}), {ext: 'sdf'});
        component.addRepresentation('ball+stick', {multipleBond: true, colorScheme: color_scheme});
        component.autoView();
        output.value = result.log || 'Molecule generated.';
    } catch (error) {
        stage.removeAllComponents();
        output.value = error.message || 'The molecule could not be displayed.';
    } finally {
        buttons.forEach(button => { button.disabled = false; button.classList.remove('grayed'); });
    }
}

function start_smiles_viewer() {
    try {
        set_stage();
        go(document.getElementById('smiles').value);
    } catch (error) {
        document.getElementById('result').value = 'Viewer startup failed: ' + error.message;
    }
}
