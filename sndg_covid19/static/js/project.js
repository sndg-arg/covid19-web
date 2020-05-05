/* Project specific Javascript goes here. */
import msa from 'msa';
import Phylocanvas from 'phylocanvas';
import Sequence from 'sequence-viewer';
import $ from 'jquery';
import 'bootstrap/dist/css/bootstrap.min.css';
import FeatureViewer from 'feature-viewer';
import * as NGL from 'ngl';

import dt from  'datatables.net';
import 'datatables.net-dt/css/jquery.dataTables.min.css';
window.$.DataTable = dt;
window.NGL = NGL;


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

window.Sequence = Sequence;
window.FeatureViewer = FeatureViewer;
window.$ = $;
