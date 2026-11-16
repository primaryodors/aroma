(function (root) {
    'use strict';
    const background = '#222222';
    const colors = {text: '#e6edf3', hbond: '#94d99a', ionic: '#94d99a',
        polarpi: '#c0d99a', mcoord: '#ffd580', pi: '#dcb3ee', covalent: '#e6edf3', other: '#aab4c0'};
    const residueColors = {A: '#dedede', I: '#dedede', L: '#dedede', V: '#dedede', P: '#dedede',
        M: '#e4dc88', G: '#b7c0cb', C: '#ffe066', F: '#e1adeb', W: '#e1adeb', Y: '#f5a9f5',
        S: '#6ee7ef', T: '#6ee7ef', N: '#6ee7ef', Q: '#6ee7ef', D: '#ff9f9f', E: '#ff9f9f',
        K: '#86baff', R: '#86baff', H: '#c4b3ff'};
    function rgb(color) {
        if (/^#[0-9a-f]{3}$/i.test(color)) color = '#' + [...color.slice(1)].map(c => c+c).join('');
        if (/^#[0-9a-f]{6}$/i.test(color)) return [1,3,5].map(i => parseInt(color.slice(i,i+2),16));
        const match = color.match(/^rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/);
        return match ? match.slice(1).map(Number) : null;
    }
    function luminance(color) {
        const channels = rgb(color).map(c => c/255).map(c => c <= .04045 ? c/12.92 : ((c+.055)/1.055)**2.4);
        return channels[0]*.2126 + channels[1]*.7152 + channels[2]*.0722;
    }
    function contrast(a,b) { const x=luminance(a), y=luminance(b); return (Math.max(x,y)+.05)/(Math.min(x,y)+.05); }
    function readable(color) {
        const value = rgb(color);
        if (!value) return color === 'none' ? color : colors.text;
        if (Math.max(...value)<30) return colors.text;
        for (let step=0; step<=20; step++) {
            const mixed = '#' + value.map(c => Math.round(c+(255-c)*step/20).toString(16).padStart(2,'0')).join('');
            if (contrast(mixed,background)>=4.5) return mixed;
        }
        return '#ffffff';
    }
    function parseMolfile(sdf) {
        const lines=sdf.replace(/\r/g,'').split('\n');
        if (!sdf.trim()) throw new Error('The cached ligand structure is unavailable.');
        if (!lines[3] || !lines[3].includes('V2000')) throw new Error('Unsupported ligand structure format.');
        const count=Number(lines[3].slice(0,3)), bonds=Number(lines[3].slice(3,6));
        if (!count || count>999) throw new Error('Invalid ligand atom inventory.');
        const atoms=[];
        for (let i=0;i<count;i++) {
            const line=lines[4+i].padEnd(69), element=line.slice(31,34).trim();
            atoms.push({name: element+(i+1), element, neighbors: []});
            lines[4+i]=line.slice(0,60)+String(i+1).padStart(3)+line.slice(63);
        }
        const edges=[];
        for (let i=0;i<bonds;i++) {
            const line=lines[4+count+i], fields=line.slice(6).trim().split(/\s+/), a=Number(line.slice(0,3))-1, b=Number(line.slice(3,6))-1;
            if (!atoms[a] || !atoms[b]) throw new Error('Invalid ligand bond.');
            const order=Number(fields[0])===1.5?4:Number(fields[0]);
            if (!Number.isInteger(order) || order<1 || order>4) throw new Error('Invalid ligand bond order.');
            const stereo=Number(fields[1]||0);
            lines[4+count+i]=[a+1,b+1,order,stereo,0,0,0].map(n=>String(n).padStart(3)).join('');
            atoms[a].neighbors.push(b); atoms[b].neighbors.push(a); edges.push([a,b]);
        }
        return {molfile: lines.join('\n'), atoms, edges};
    }
    const median = values => { const list=values.slice().sort((a,b)=>a-b); return list.length ? list[Math.floor(list.length/2)] : 1; };
    function align(points, projected, edges) {
        const names=Object.keys(points).filter(name => projected[name] && !/^H\d+$/.test(name));
        if (!names.length) throw new Error('Ligand atoms do not match this dock.');
        const center = source => names.reduce((sum,name) => [sum[0]+source[name][0]/names.length, sum[1]+source[name][1]/names.length], [0,0]);
        const from=center(points), to=center(projected);
        let dot=0, cross=0;
        for (const name of names) {
            const x=points[name][0]-from[0], y=points[name][1]-from[1];
            const u=projected[name][0]-to[0], v=projected[name][1]-to[1];
            dot+=x*u+y*v; cross+=x*v-y*u;
        }
        const clean=[], docked=[];
        for (const [a,b] of edges) if (points[a] && points[b] && projected[a] && projected[b] && !/^H\d+$/.test(a) && !/^H\d+$/.test(b)) {
            clean.push(Math.hypot(points[a][0]-points[b][0],points[a][1]-points[b][1]));
            const distance=Math.hypot(projected[a][0]-projected[b][0],projected[a][1]-projected[b][1]);
            if (distance>1) docked.push(distance);
        }
        const scale=Math.max(.5,Math.min(2,median(docked)/median(clean)));
        const angle=Math.atan2(cross,dot), a=scale*Math.cos(angle), b=scale*Math.sin(angle);
        return [a,b,-b,a,to[0]-a*from[0]+b*from[1],to[1]-b*from[0]-a*from[1]];
    }
    const transform = (p,m) => [m[0]*p[0]+m[2]*p[1]+m[4], m[1]*p[0]+m[3]*p[1]+m[5]];
    function svgAtomPoints(svg,names) {
        const points={};
        for (const tag of svg.match(/<circle\b[^>]*>/g)||[]) {
            const index=tag.match(/:Atom:(\d+)"/), x=tag.match(/\bcx="([^"]+)"/), y=tag.match(/\bcy="([^"]+)"/);
            if (index && x && y) points[names[Number(index[1])]]=[Number(x[1]),Number(y[1])];
        }
        return points;
    }
    function xmlSafeSVG(svg) {
        // OpenChemLib 5 emits relative stereochemistry labels such as &1 without XML escaping.
        return svg.replace(/&(?!(?:amp|lt|gt|quot|apos|#\d+|#x[\da-fA-F]+);)/g,'&amp;')
            // Bond geometry already conveys E/Z; retain the atom R/S descriptors.
            .replace(/<text\b[^>]*>\s*[EZ]\s*<\/text>/g,'');
    }
    function cachedFallback(smiles,positions,OCL) {
        const molecule=OCL.Molecule.fromSmiles(smiles);
        const radii={C:.76,N:.71,O:.66,S:1.05,P:1.07,F:.57,Cl:1.02,Br:1.2,I:1.39};
        const names=Object.keys(positions).filter(name=>!/^H\d+$/.test(name)).sort((a,b)=>Number(a.match(/\d+$/)[0])-Number(b.match(/\d+$/)[0]));
        if(names.length!==molecule.getAllAtoms()) throw new Error('Catalog and dock atom inventories differ.');
        for(let i=0;i<names.length;i++) {
            if(names[i]!==molecule.getAtomLabel(i)+(i+1)) throw new Error('Catalog atom order differs from this dock.');
            for(let j=0;j<i;j++) {
                const a=molecule.getAtomLabel(i), b=molecule.getAtomLabel(j);
                if(!radii[a] || !radii[b]) throw new Error('Catalog atom mapping cannot be verified.');
                const distance=Math.hypot(...positions[names[i]].map((v,k)=>v-positions[names[j]][k]));
                if((distance<radii[a]+radii[b]+.35)!==(molecule.getBond(i,j)>=0))
                    throw new Error('Catalog bonds do not match this dock.');
            }
        }
        for(const name of Object.keys(positions).filter(name=>/^H\d+$/.test(name)).sort((a,b)=>Number(a.slice(1))-Number(b.slice(1)))) {
            const h=molecule.addAtom(1);
            if(name!=='H'+(h+1)) throw new Error('Dock hydrogen inventory differs from the catalog.');
            const distances=names.map(atom=>Math.hypot(...positions[atom].map((v,k)=>v-positions[name][k])));
            const parent=distances.indexOf(Math.min(...distances));
            if(distances[parent]>1.5) throw new Error('Dock hydrogen parent cannot be verified.');
            molecule.addBond(parent,h,1);
        }
        return molecule.toMolfile();
    }
    function depiction(sdf, projected, interactions, OCL, id) {
        const data=parseMolfile(sdf), molecule=OCL.Molecule.fromMolfile(data.molfile);
        if (molecule.getAllBonds()!==data.edges.length) throw new Error('Ligand bonds were not decoded correctly.');
        const linkAtoms={};
        for (let i=0;i<data.atoms.length;i++) {
            const atom=data.atoms[i];
            if (atom.element!=='O' || projected[atom.name] || atom.neighbors.length!==1) continue;
            const carbon=atom.neighbors[0], bond=molecule.getBond(i,carbon);
            const contact=interactions.find(item=>item.type==='covalent' && item.atom===data.atoms[carbon].name
                && /^N/.test(item.proteinAtom||''));
            if (!contact || data.atoms[carbon].element!=='C' || molecule.getBondOrder(bond)!==2) continue;
            atom.name='protein:'+contact.residue+':'+contact.proteinAtom;
            atom.external=true;
            molecule.setAtomicNo(i,7);
            // The remaining nitrogen valence connects to the protein outside this fragment.
            molecule.setAtomAbnormalValence(i,2);
            linkAtoms[contact.residue+':'+contact.atom]=atom.name;
        }
        molecule.inventCoordinates();
        molecule.ensureHelperArrays(OCL.Molecule.cHelperSymmetrySimple);
        for (let b=0;b<molecule.getAllBonds();b++) {
            if (molecule.getBondType(b)!==OCL.Molecule.cBondTypeCross) continue;
            for (const side of [0,1]) {
                const atom=molecule.getBondAtom(side,b), other=molecule.getBondAtom(1-side,b), branches=[];
                for(let j=0;j<molecule.getConnAtoms(atom);j++) {
                    const neighbor=molecule.getConnAtom(atom,j);
                    if(neighbor!==other) branches.push(molecule.getSymmetryRank(neighbor));
                }
                if (molecule.getAtomicNo(atom)===6 && ((branches.length===2 && branches[0]===branches[1])
                    || molecule.getImplicitHydrogens(atom)===2)) molecule.setBondType(b,OCL.Molecule.cBondTypeDouble);
            }
        }
        const indexByOriginal = () => new Map(Array.from({length:molecule.getAllAtoms()},(_,i)=>[molecule.getAtomMapNo(i)-1,i]));
        let indices=indexByOriginal();
        for (const name of new Set(interactions.map(item=>item.atom))) {
            const original=data.atoms.findIndex(atom=>atom.name===name);
            if (original<0 || data.atoms[original].element!=='H' || indices.has(original)) continue;
            const parent=data.atoms[original].neighbors.find(i=>data.atoms[i].element!=='H');
            const current=indices.get(parent);
            if (current===undefined) throw new Error('Interacting hydrogen has no ligand parent.');
            const origin=[molecule.getAtomX(current),molecule.getAtomY(current)];
            const angles=[];
            for (let i=0;i<molecule.getAllAtoms();i++) if (i!==current && molecule.getBond(current,i)>=0)
                angles.push(Math.atan2(molecule.getAtomY(i)-origin[1],molecule.getAtomX(i)-origin[0]));
            angles.sort((a,b)=>a-b);
            let gap=-1, angle=0;
            for (let i=0;i<angles.length;i++) {
                const next=i+1<angles.length?angles[i+1]:angles[0]+Math.PI*2;
                if (next-angles[i]>gap) {gap=next-angles[i]; angle=(next+angles[i])/2;}
            }
            const h=molecule.addAtom(1); molecule.addBond(current,h,1);
            molecule.setAtomX(h,origin[0]+.8*Math.cos(angle)); molecule.setAtomY(h,origin[1]+.8*Math.sin(angle));
            molecule.setAtomMapNo(h,original+1,false);
            indices=indexByOriginal();
        }
        molecule.setHydrogenProtection(true); molecule.ensureHelperArrays(1);
        const names=Array.from({length:molecule.getAllAtoms()},(_,i)=>data.atoms[molecule.getAtomMapNo(i)-1]?.name);
        if (names.some(name=>!name)) throw new Error('Ligand atom mapping was lost.');
        for (let i=0;i<molecule.getAllAtoms();i++) molecule.setAtomMapNo(i,0,false);
        const options={fontWeight:600,suppressCIPParity:false,suppressChiralText:true,suppressESR:true};
        let svg=xmlSafeSVG(molecule.toSVG(300,300,id,options)), points=svgAtomPoints(svg,names);
        for (const atom of data.atoms) if (atom.element!=='H' && !atom.external && (!points[atom.name] || !projected[atom.name]))
            throw new Error('Ligand atom names do not match this dock.');
        const edges=data.edges.map(([a,b])=>[data.atoms[a].name,data.atoms[b].name]);
        const initial=align(points,projected,edges), angle=Math.atan2(initial[1],initial[0]);
        // Rotate molecular coordinates before depiction so atom labels stay upright.
        for (let i=0;i<molecule.getAllAtoms();i++) {
            const x=molecule.getAtomX(i), y=molecule.getAtomY(i);
            molecule.setAtomX(i,x*Math.cos(angle)-y*Math.sin(angle));
            molecule.setAtomY(i,x*Math.sin(angle)+y*Math.cos(angle));
        }
        svg=xmlSafeSVG(molecule.toSVG(300,300,id,options)); points=svgAtomPoints(svg,names);
        const aligned=align(points,projected,edges), scale=Math.hypot(aligned[0],aligned[1]);
        const matched=Object.keys(points).filter(name=>projected[name] && !/^H\d+$/.test(name));
        const offset=matched.reduce((sum,name)=>[sum[0]+(projected[name][0]-scale*points[name][0])/matched.length,
            sum[1]+(projected[name][1]-scale*points[name][1])/matched.length],[0,0]);
        const matrix=[scale,0,0,scale,...offset];
        const anchors=Object.fromEntries(Object.entries(points).map(([name,point])=>[name,transform(point,matrix)]));
        const ions=[];
        for (const [name,point] of Object.entries(projected))
            if (!anchors[name] && /^(Cu|Zn|Fe|Co|Ni|Mn|Mg|Ca|Na|K|Ag|Au|Cd|Hg|Pt|Pd)\d+$/.test(name)) {
                anchors[name]=point; ions.push({name,point,label:name.replace(/\d+/g,'')});
            }
        for (const item of interactions) if (!anchors[item.atom]) throw new Error('Interaction atom is absent from the ligand depiction: '+item.atom);
        return {svg,matrix,anchors,names,ions,linkAtoms};
    }
    function render(host, data) {
        const ns='http://www.w3.org/2000/svg';
        const node=(tag,attributes={},text) => {const result=document.createElementNS(ns,tag); for(const [key,value] of Object.entries(attributes)) result.setAttribute(key,value); if(text!==undefined) result.textContent=text; return result;};
        try {
            data={...data,residues:Object.fromEntries(Object.entries(data.residues).map(([key,value])=>[key,{...value}]))};
            const source=data.sdf || cachedFallback(data.smiles,data.positions,root.OCL);
            const id=host.id+'-ligand', drawn=depiction(source,data.projected,data.interactions,root.OCL,id);
            const svg=node('svg',{xmlns:ns,role:'img','aria-labelledby':host.id+'-title '+host.id+'-description'});
            svg.append(node('title',{id:host.id+'-title'},data.title));
            svg.append(node('desc',{id:host.id+'-description'},'Chemical ligand depiction with receptor residues at their dock diagram positions. '+data.interactions.map(i=>i.residue+' to '+i.atom+': '+i.type).join('; ')));
            const drawInteractions=() => { for (const item of data.interactions) {
                const residue=data.residues[item.residue]; if (!residue) continue;
                const start=drawn.anchors[drawn.linkAtoms[item.residue+':'+item.atom]||item.atom], distance=Math.hypot(residue.x-start[0],residue.y-start[1]);
                if (distance<1) continue;
                const end=[residue.x+(start[0]-residue.x)*25/distance,residue.y+(start[1]-residue.y)*25/distance];
                const patterns={hbond:'6 3',ionic:'2 3',pi:'7 3 2 3',polarpi:'5 3',mcoord:'5 3',covalent:''};
                const line=node('line',{x1:start[0],y1:start[1],x2:end[0],y2:end[1],stroke:colors[item.type]||colors.other,
                    'stroke-width':1.8,'stroke-dasharray':patterns[item.type]??'2 5',class:'show'+item.residue.replace('.','x'),
                    'data-ligand-atom':item.atom,'data-residue':item.residue});
                const linked=drawn.linkAtoms[item.residue+':'+item.atom];
                if(linked) line.setAttribute('data-protein-atom',item.proteinAtom);
                line.append(node('title',{},linked ? item.atom+'=N – '+item.residue+':'+item.proteinAtom+' (Schiff linkage)'
                    : item.residue+' – '+item.atom+': '+item.type)); svg.insertBefore(line,svg.querySelector('.dock-ligand'));
            }};
            const group=node('g',{transform:'matrix('+drawn.matrix.join(' ')+')',class:'dock-ligand'});
            const parsed=new DOMParser().parseFromString(drawn.svg,'image/svg+xml').documentElement;
            if (parsed.localName!=='svg') throw new Error('The ligand SVG could not be read.');
            const allowed=new Set(['g','line','path','polygon','polyline','circle','ellipse','text','tspan','rect']);
            for (const element of [...parsed.querySelectorAll('*')]) {
                if (!allowed.has(element.localName)) { element.remove(); continue; }
                for (const attribute of [...element.attributes]) if (/^on/i.test(attribute.name)||/href/i.test(attribute.name)) element.removeAttribute(attribute.name);
                for (const property of ['fill','stroke']) {
                    const color=element.style[property] || element.getAttribute(property);
                    if (color) { element.style[property]=readable(color); element.removeAttribute(property); }
                }
                if (element.localName==='line' && !element.classList.contains('event')) element.style.strokeWidth='1.5';
                if (element.localName==='text') { element.style.paintOrder='stroke'; element.style.stroke=background; element.style.strokeWidth='3px'; element.style.strokeLinejoin='round'; }
                const atom=element.id.match(/:Atom:(\d+)$/);
                if (atom) {const name=drawn.names[Number(atom[1])]; element.setAttribute('data-atom-name',name); element.append(node('title',{},name));}
            }
            while(parsed.firstChild) group.append(parsed.firstChild); svg.append(group);
            for (const ion of drawn.ions) {
                const label=node('text',{x:ion.point[0],y:ion.point[1],fill:colors.mcoord,'font-size':16,'text-anchor':'middle'},ion.label);
                label.append(node('title',{},ion.name)); svg.append(label);
            }
            host.replaceChildren(svg);
            const obstacles=[];
            for (const element of group.querySelectorAll('line,text,path,polygon,polyline')) {
                if(element.classList.contains('event')) continue;
                const box=element.getBBox(), m=drawn.matrix;
                const a=transform([box.x,box.y],m), b=transform([box.x+box.width,box.y+box.height],m);
                obstacles.push([a[0]-5,a[1]-5,b[0]+5,b[1]+5]);
            }
            const placed=[];
            for (const residue of Object.values(data.residues)) {
                const original=[residue.x,residue.y];
                const clear=(x,y)=>obstacles.every(b=>Math.hypot(x-Math.max(b[0],Math.min(x,b[2])),y-Math.max(b[1],Math.min(y,b[3])))>33)
                    && placed.every(p=>Math.hypot(x-p[0],y-p[1])>56);
                if(!clear(...original)) {
                    search: for(let radius=4;radius<2000;radius+=4) for(let step=0;step<64;step++) {
                        const angle=step*Math.PI/32, x=original[0]+radius*Math.cos(angle), y=original[1]+radius*Math.sin(angle);
                        if(clear(x,y)) {residue.x=x;residue.y=y;break search;}
                    }
                }
                placed.push([residue.x,residue.y]);
            }
            drawInteractions();
            for (const [bw,residue] of Object.entries(data.residues)) {
                const color=residueColors[residue.letter]||colors.text, g=node('g',{class:'show'+bw.replace('.','x')});
                g.append(node('circle',{cx:residue.x,cy:residue.y,r:25,fill:background,stroke:color,'stroke-width':1.5}));
                g.append(node('text',{x:residue.x,y:residue.y-3,fill:color,'text-anchor':'middle','font-size':12},bw));
                g.append(node('text',{x:residue.x,y:residue.y+11,fill:color,'text-anchor':'middle','font-size':12},residue.letter+residue.number)); svg.append(g);
            }
            const points=[...Object.values(drawn.anchors),...Object.values(data.residues).map(r=>[r.x,r.y])];
            const x=Math.min(0,...points.map(p=>p[0]-36)), y=Math.min(0,...points.map(p=>p[1]-36));
            const width=Math.max(data.width,...points.map(p=>p[0]+36))-x, height=Math.max(data.height,...points.map(p=>p[1]+36))-y;
            svg.setAttribute('viewBox',[x,y,width,height].join(' ')); svg.setAttribute('width',width); svg.setAttribute('height',height);
            svg.insertBefore(node('rect',{x,y,width,height,fill:background}),svg.children[2]||null);
            host.replaceChildren(svg);
            const legend=document.createElement('p'); legend.className='dock-diagram-legend';
            legend.textContent='Interactions: long dashes — hydrogen bonds; short dashes — ionic contacts; dash-dot — π contacts; medium dashes — metal/polar–π contacts; fine dots — other contacts; solid — covalent contacts. Hover for atom names and interaction types.';
            host.append(legend);
        } catch (error) {
            host.textContent='2D ligand diagram unavailable: '+error.message; host.setAttribute('role','status');
        }
    }
    const api={background,colors,residueColors,contrast,readable,parseMolfile,align,transform,cachedFallback,depiction,render};
    if(typeof module!=='undefined' && module.exports) module.exports=api;
    root.DockDiagram=api;
})(typeof window==='undefined'?globalThis:window);
