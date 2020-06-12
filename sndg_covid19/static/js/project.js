/* Project specific Javascript goes here. */
import msa from 'msa';
import Phylocanvas from 'phylocanvas';
import Sequence from 'sequence-viewer';
import $ from 'jquery';
import 'bootstrap/dist/css/bootstrap.min.css';
import FeatureViewer from 'feature-viewer';
import * as NGL from 'ngl';

import 'pivottable/dist/pivot.css'
import 'pivottable';


import dt from  'datatables.net';
import 'datatables.net-dt/css/jquery.dataTables.min.css';
window.$.DataTable = dt;
window.NGL = NGL;

window.msa = msa;


// import 'bootstrap/dist/css/bootstrap.css';

class SeqDiv extends HTMLElement {
    constructor() {
        // Always call super first in constructor
        super();

        // write element functionality in here


        const seq = new Sequence( "RFQAEGSLKK");
        const properties =  {
            'showLineNumbers': true,
            'wrapAminoAcids': true,
            'charsPerLine':  100,
            'toolbar': false,
            'search': false,
            'title': "a seq",
            //'sequenceMaxHeight': "30px",
            'badge': false
        };
        seq.render( '#ajujujaja' , properties);
    }
}



window.customElements.define("sequence-viewer", SeqDiv);

// import  '@gmod/jbrowse/dist/main.bundle';

// window.JBrowse = JBrowse;

window.Sequence = Sequence;
window.FeatureViewer = FeatureViewer;
window.$ = $;


