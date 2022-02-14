var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __decorate = (this && this.__decorate) || function (decorators, target, key, desc) {
    var c = arguments.length, r = c < 3 ? target : desc === null ? desc = Object.getOwnPropertyDescriptor(target, key) : desc, d;
    if (typeof Reflect === "object" && typeof Reflect.decorate === "function") r = Reflect.decorate(decorators, target, key, desc);
    else for (var i = decorators.length - 1; i >= 0; i--) if (d = decorators[i]) r = (c < 3 ? d(r) : c > 3 ? d(target, key, r) : d(target, key)) || r;
    return c > 3 && r && Object.defineProperty(target, key, r), r;
};
import { dialogAddons, SortByDefault, toolbar } from "lineupjs";
import MapColumn from "lineupjs";
import { isMissingValue } from "lineupjs";
/**
 * a numbersMap column with optional alignment
 */
var NumbersMapColumn = /** @class */ (function (_super) {
    __extends(NumbersMapColumn, _super);
    function NumbersMapColumn(id, desc) {
        var _this = _super.call(this, id, desc) || this;
        _this.dataLength = desc.dataLength;
        return _this;
    }
    NumbersMapColumn.prototype.on = function (type, listener) {
        return _super.prototype.on.call(this, type, listener);
    };
    NumbersMapColumn.prototype.getValue = function (row) {
        var r = this.getMapValue(row);
        return r.every(function (d) { return d.value === ""; }) ? null : r;
    };
    NumbersMapColumn.prototype.getMapValue = function (row) {
        return _super.prototype.getMap.call(this, row).map(function (_a) {
            var key = _a.key, value = _a.value;
            return ({
                key: key,
                value: isMissingValue(value) ? "" : String(value),
            });
        });
    };
    NumbersMapColumn = __decorate([
        toolbar("rename", "filterNumber", "colorMapped", "editMapping"),
        dialogAddons("sort", "sortNumbers"),
        SortByDefault("descending")
        //@ts-ignore
    ], NumbersMapColumn);
    return NumbersMapColumn;
}(MapColumn));
export default NumbersMapColumn;
//# sourceMappingURL=NumbersMapColumn.js.map