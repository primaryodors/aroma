const assert = require('node:assert/strict');
const fs = require('node:fs');
const diagram = require('../www/js/dock_diagram.js');
assert.ok(Math.abs(diagram.contrast('#000000','#ffffff')-21)<1e-12);
assert.equal(diagram.contrast('#222222','#222222'),1);

for (const color of [...Object.values(diagram.colors), ...Object.values(diagram.residueColors)])
    assert.ok(diagram.contrast(color,diagram.background)>=4.5, color);
for (const backdrop of ['#222','#245','#210','#444','#111','#226','#520','#232','#313','#293541'])
    for (const color of [...Object.values(diagram.residueColors),'#9acbff','#e6edf3'])
        assert.ok(diagram.contrast(color,backdrop)>=4.5, color+' on '+backdrop);
for (const color of ['rgb(255,13,13)','rgb(0,0,255)','rgb(160,0,0)','#000','#c60','#609'])
    assert.ok(diagram.contrast(diagram.readable(color),diagram.background)>=4.5);
const stylesheet=fs.readFileSync('www/assets/dock_diagram.css','utf8');
for (const match of stylesheet.matchAll(/(?<![-\w])color:\s*(#[0-9a-f]{6})/gi))
    assert.ok(diagram.contrast(match[1],'#444444')>=4.5,match[1]);

const points={C1:[0,0],C2:[1,0],O3:[1,1]}, projected={C1:[10,20],C2:[10,21.5],O3:[8.5,21.5]};
const matrix=diagram.align(points,projected,[['C1','C2'],['C2','O3']]);
assert.ok(matrix[0]*matrix[3]-matrix[1]*matrix[2]>0);
for (const name of Object.keys(points)) {
    const actual=diagram.transform(points[name],matrix);
    assert.ok(Math.hypot(actual[0]-projected[name][0],actual[1]-projected[name][1])<1e-8);
}

if (process.argv[2]) {
    const OCL=require(process.argv[2]);
    const sources=['sdf/phenethyl_alcohol.sdf','testdata/benzyl_salicylate_contact.sdf','testdata/hydroxy_acid_contact.sdf',
        'sdf/geraniol.sdf','sdf/cinnamaldehyde.sdf','sdf/peony_acetonitrile.sdf','sdf/farnesol.sdf','sdf/octanal.sdf','sdf/rosyfolia.sdf'];
    for (const filename of sources) {
        const sdf=fs.readFileSync(filename,'utf8'), parsed=diagram.parseMolfile(sdf), rows=sdf.split('\n');
        const projected=Object.fromEntries(parsed.atoms.map((atom,i)=>[atom.name,
            [300+Number(rows[4+i].slice(0,10))*25,200+Number(rows[4+i].slice(10,20))*25]]));
        const heavy=parsed.atoms.filter(atom=>atom.element!=='H');
        const drawn=diagram.depiction(sdf,projected,heavy.map(atom=>({atom:atom.name})),OCL,'test');
        for (const atom of heavy) assert.ok(drawn.anchors[atom.name],atom.name);
        assert.ok(!/&(?!(?:amp|lt|gt|quot|apos|#\d+|#x[\da-fA-F]+);)/.test(drawn.svg),filename);
        if(filename.includes('rosyfolia')) {
            assert.ok(!drawn.svg.includes('&amp;1'),'Relative stereo group annotations are hidden.');
            assert.equal((drawn.svg.match(/>R<\/text>/g)||[]).length,2);
            assert.equal((drawn.svg.match(/>S<\/text>/g)||[]).length,1);
            assert.ok(!drawn.svg.includes('rgb(0,128,0)'),'Stereo bonds use the normal skeletal bond color.');
            assert.ok(drawn.svg.includes('<polygon'),'Stereochemical wedges remain present.');
        }
        assert.ok(!/<text\b[^>]*>\s*[EZ]\s*<\/text>/.test(drawn.svg),filename);
        assert.throws(()=>diagram.depiction(sdf,projected,[{atom:'C999'}],OCL,'missing'),/absent/);
        if (filename.includes('octanal')) {
            const reacted={...projected}; delete reacted.O9;
            assert.throws(()=>diagram.depiction(sdf,reacted,[],OCL,'unexplained'),/do not match/);
            assert.throws(()=>diagram.depiction(sdf,reacted,[{atom:'C8',type:'covalent',proteinAtom:'SG',residue:'4.60'}],OCL,'other'),/do not match/);
            const adduct=diagram.depiction(sdf,reacted,[{atom:'C8',type:'covalent',proteinAtom:'NZ',residue:'4.60'}],OCL,'schiff');
            assert.equal(adduct.anchors.O9,undefined);
            const nitrogen=adduct.linkAtoms['4.60:C8'];
            assert.equal(nitrogen,'protein:4.60:NZ');
            assert.ok(adduct.anchors[nitrogen]);
            assert.ok(adduct.svg.includes('>N</text>'));
            const c=adduct.names.indexOf('C8'), n=adduct.names.indexOf(nitrogen);
            assert.ok(adduct.svg.includes('Bond:'+c+'-'+n) || adduct.svg.includes('Bond:'+n+'-'+c));
            assert.equal((adduct.svg.match(/<line x1=/g)||[]).length,9,'Seven chain bonds and two imine strokes.');
        }
        if (filename.includes('phenethyl')) {
            const positions=Object.fromEntries(parsed.atoms.map((atom,i)=>[atom.name,
                [0,10,20].map(offset=>Number(rows[4+i].slice(offset,offset+10)))]));
            const fallback=diagram.parseMolfile(diagram.cachedFallback('c1ccc(cc1)CCO',positions,OCL));
            assert.deepEqual(fallback.atoms.map(a=>a.name),parsed.atoms.map(a=>a.name));
            assert.deepEqual(fallback.edges.map(e=>e.slice().sort().join('-')).sort(),parsed.edges.map(e=>e.slice().sort().join('-')).sort());
            assert.throws(()=>diagram.cachedFallback('c1ccc(cc1)CCO',{...positions,O9:[100,100,100]},OCL),/bonds do not match/);
            const polar=diagram.depiction(sdf,projected,[{atom:'H19'}],OCL,'polar');
            assert.ok(polar.anchors.H19);
            const h=polar.names.indexOf('H19'), o=polar.names.indexOf('O9');
            assert.ok(polar.svg.includes('Bond:'+o+'-'+h) || polar.svg.includes('Bond:'+h+'-'+o));
            const lengths=Array.from({length:6},(_,i)=>Math.hypot(
                polar.anchors['C'+(i+1)][0]-polar.anchors['C'+((i+1)%6+1)][0],
                polar.anchors['C'+(i+1)][1]-polar.anchors['C'+((i+1)%6+1)][1]));
            assert.ok(Math.max(...lengths)/Math.min(...lengths)<1.1);
            const metal=diagram.depiction(sdf,{...projected,Cu1:[20,30]},[{atom:'Cu1'}],OCL,'metal');
            assert.deepEqual(metal.anchors.Cu1,[20,30]);
        }
    }
    console.log('Real OpenChemLib atom identity, hydrogen endpoints, regular rings, and metal anchors pass.');
}
console.log('Contrast ratios, residue palettes, and rotation without reflection pass.');
