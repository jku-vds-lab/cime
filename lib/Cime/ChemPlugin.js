var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
import * as React from 'react';
import { DatasetType, PSEPlugin } from 'projection-space-explorer';
import { ChemLegendParent } from './ChemDetail/ChemDetail';
var ChemPlugin = /** @class */ (function (_super) {
    __extends(ChemPlugin, _super);
    function ChemPlugin() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.type = DatasetType.Chem;
        return _this;
    }
    // @ts-ignore
    ChemPlugin.prototype.createFingerprint = function (dataset, vectors, scale, aggregate) {
        return React.createElement(ChemLegendParent, { selection: vectors, aggregate: aggregate });
    };
    return ChemPlugin;
}(PSEPlugin));
export { ChemPlugin };
//# sourceMappingURL=ChemPlugin.js.map