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
var __decorate = (this && this.__decorate) || function (decorators, target, key, desc) {
    var c = arguments.length, r = c < 3 ? target : desc === null ? desc = Object.getOwnPropertyDescriptor(target, key) : desc, d;
    if (typeof Reflect === "object" && typeof Reflect.decorate === "function") r = Reflect.decorate(decorators, target, key, desc);
    else for (var i = decorators.length - 1; i >= 0; i--) if (d = decorators[i]) r = (c < 3 ? d(r) : c > 3 ? d(target, key, r) : d(target, key)) || r;
    return c > 3 && r && Object.defineProperty(target, key, r), r;
};
import { format } from "d3";
import { Column, dialogAddons, EAdvancedSortMethod, ECompareValueType, MapColumn, NumberColumn, ScaleMappingFunction, SortByDefault, toolbar, } from "lineupjs";
import { DEFAULT_FORMATTER, isDummyNumberFilter, noNumberFilter, restoreMapping, restoreNumberFilter, } from "./helper_methods";
var TestColumn = /** @class */ (function (_super) {
    __extends(TestColumn, _super);
    function TestColumn(id, desc, factory) {
        var _this = _super.call(this, id, desc) || this;
        _this.numberFormat = DEFAULT_FORMATTER;
        /**
         * currently active filter
         * @type {{min: number, max: number}}
         * @private
         */
        _this.currentFilter = noNumberFilter();
        _this.min = 0;
        _this.max = 1;
        // this.mapping = restoreMapping(desc, factory); // TODO: check, if desc.range and desc.domain can be infered
        _this.mapping = new ScaleMappingFunction([desc["min"], desc["max"]], "linear", [0, 1]);
        _this.original = _this.mapping.clone();
        _this.sort = desc.sort || EAdvancedSortMethod.median;
        _this.colorMapping = factory.colorMappingFunction(desc.colorMapping || desc.color);
        if (desc.numberFormat) {
            _this.numberFormat = format(desc.numberFormat);
        }
        //TODO: infer min and max if it is not given
        _this.min = desc["min"];
        _this.max = desc["max"];
        return _this;
    }
    TestColumn_1 = TestColumn;
    TestColumn.prototype.getMin = function () {
        return this.min;
    };
    TestColumn.prototype.getMax = function () {
        return this.max;
    };
    TestColumn.prototype.getNumberFormat = function () {
        return this.numberFormat;
    };
    // https://stackoverflow.com/questions/45309447/calculating-median-javascript
    TestColumn.prototype.get_quartile = function (values, q) {
        if (q === void 0) { q = 0.5; }
        // 1. quartile: q=0.25 | median: q=0.5 | 3. quartile: q=0.75
        if (values.length === 0)
            return 0;
        values.sort(function (a, b) {
            return a - b;
        });
        var half = Math.floor(values.length * q);
        if (values.length % 2)
            return values[half];
        return (values[half - 1] + values[half]) / 2.0;
    };
    // https://www.sitepoint.com/community/t/calculating-the-average-mean/7302/2
    TestColumn.prototype.mean = function (numbers) {
        var total = 0, i;
        for (i = 0; i < numbers.length; i += 1) {
            total += numbers[i];
        }
        return total / numbers.length;
    };
    TestColumn.prototype.get_advanced_value = function (method, value_list) {
        switch (method) {
            case EAdvancedSortMethod.min:
                return Math.min.apply(Math, value_list);
            case EAdvancedSortMethod.max:
                return Math.max.apply(Math, value_list);
            case EAdvancedSortMethod.mean:
                return this.mean(value_list);
            case EAdvancedSortMethod.median:
                return this.get_quartile(value_list);
            case EAdvancedSortMethod.q1:
                return this.get_quartile(value_list, 1);
            case EAdvancedSortMethod.q3:
                return this.get_quartile(value_list, 3);
            default:
                return this.get_quartile(value_list);
        }
    };
    TestColumn.prototype.toCompareValue = function (row) {
        var data = this.getValue(row);
        var value_list = data[0]["value"];
        var method = this.getSortMethod();
        return this.get_advanced_value(method, value_list);
    };
    TestColumn.prototype.toCompareValueType = function () {
        return ECompareValueType.FLOAT;
    };
    TestColumn.prototype.getBoxPlotDataFromValueList = function (data) {
        return {
            mean: this.get_advanced_value(EAdvancedSortMethod.mean, data),
            missing: 0,
            count: data.length,
            max: this.get_advanced_value(EAdvancedSortMethod.max, data),
            min: this.get_advanced_value(EAdvancedSortMethod.min, data),
            median: this.get_advanced_value(EAdvancedSortMethod.median, data),
            q1: this.get_advanced_value(EAdvancedSortMethod.q1, data),
            q3: this.get_advanced_value(EAdvancedSortMethod.q3, data),
        };
    };
    TestColumn.prototype.getBoxPlotData = function (row) {
        console.log("getBoxPlotData");
        var data = this.getValue(row)[0]["value"];
        if (data == null) {
            return null;
        }
        return this.getBoxPlotDataFromValueList(data);
    };
    TestColumn.prototype.getRawBoxPlotData = function (row) {
        console.log("getRawBoxPlotData");
        var data = this.getRawValue(row)[0]["value"];
        if (data == null) {
            return null;
        }
        return this.getBoxPlotDataFromValueList(data);
    };
    TestColumn.prototype.getRange = function () {
        console.log("getRange");
        return this.mapping.getRange(this.numberFormat);
    };
    TestColumn.prototype.getColorMapping = function () {
        console.log("getColorMapping");
        return this.colorMapping.clone();
    };
    TestColumn.prototype.getNumber = function (row) {
        // console.log("getNumber")
        return this.mapping.apply(this.toCompareValue(row));
    };
    TestColumn.prototype.getRawNumber = function (row) {
        // console.log("getRawNumber")
        return this.toCompareValue(row);
    };
    TestColumn.prototype.iterNumber = function (row) {
        // console.log("iterNumber")
        var r = this.getValue(row);
        // return r ? r.map((d) => d.value) : [NaN];
        // return r ? r[0]["value"] : [NaN];
        return [this.get_advanced_value(EAdvancedSortMethod.median, r[0]["value"])];
    };
    TestColumn.prototype.iterRawNumber = function (row) {
        // console.log("iterRawNumber")
        var r = this.getRawValue(row);
        // return r ? r.map((d) => d.value) : [NaN];
        // return r ? r[0]["value"] : [NaN];
        return [this.get_advanced_value(EAdvancedSortMethod.median, r[0]["value"])];
    };
    TestColumn.prototype.getValue = function (row) {
        var _this = this;
        var values = this.getRawValue(row);
        if (values.length === 0) {
            //@ts-ignore
            return null;
        }
        //@ts-ignore
        return values.map(function (_a) {
            var key = _a.key, value = _a.value;
            return {
                key: key,
                value: value.length === 0
                    ? null
                    : value.map(function (val) { return _this.mapping.apply(val); }),
            };
        });
    };
    TestColumn.prototype.getRawValue = function (row) {
        var r = _super.prototype.getValue.call(this, row);
        return r == null ? [] : r;
        // const values = super.getValue(row);
        // if(values.length === 0)
        //     return null;
        // return values.map(({key, value}) => {
        //     return {key, value: value.length===0 ? null : value};
        // });
    };
    TestColumn.prototype.getExportValue = function (row, format) {
        return format === "json"
            ? this.getRawValue(row)
            : _super.prototype.getExportValue.call(this, row, format);
    };
    TestColumn.prototype.getLabels = function (row) {
        var _this = this;
        var v = this.getRawValue(row);
        return Object.keys(v).map(function (key) { return ({
            key: key,
            value: v[key].map(function (val) { return _this.numberFormat(val); }),
        }); });
    };
    TestColumn.prototype.getSortMethod = function () {
        return this.sort;
    };
    TestColumn.prototype.setSortMethod = function (sort) {
        if (this.sort === sort) {
            return;
        }
        this.fire([TestColumn_1.EVENT_SORTMETHOD_CHANGED], this.sort, (this.sort = sort));
        // sort by me if not already sorted by me
        if (!this.isSortedByMe().asc) {
            this.sortByMe();
        }
    };
    TestColumn.prototype.dump = function (toDescRef) {
        var r = _super.prototype.dump.call(this, toDescRef);
        r.sortMethod = this.getSortMethod();
        r.filter = !isDummyNumberFilter(this.currentFilter)
            ? this.currentFilter
            : null;
        r.map = this.mapping.toJSON();
        return r;
    };
    TestColumn.prototype.restore = function (dump, factory) {
        _super.prototype.restore.call(this, dump, factory);
        if (dump.sortMethod) {
            this.sort = dump.sortMethod;
        }
        if (dump.filter) {
            this.currentFilter = restoreNumberFilter(dump.filter);
        }
        if (dump.map || dump.domain) {
            this.mapping = restoreMapping(dump, factory);
        }
    };
    TestColumn.prototype.createEventList = function () {
        return _super.prototype.createEventList.call(this)
            .concat([
            TestColumn_1.EVENT_MAPPING_CHANGED,
            TestColumn_1.EVENT_SORTMETHOD_CHANGED,
            TestColumn_1.EVENT_FILTER_CHANGED,
        ]);
    };
    TestColumn.prototype.on = function (type, listener) {
        return _super.prototype.on.call(this, type, listener);
    };
    TestColumn.prototype.getOriginalMapping = function () {
        return this.original.clone();
    };
    TestColumn.prototype.getMapping = function () {
        return this.mapping.clone();
    };
    TestColumn.prototype.setMapping = function (mapping) {
        if (this.mapping.eq(mapping)) {
            return;
        }
        this.fire([
            TestColumn_1.EVENT_MAPPING_CHANGED,
            Column.EVENT_DIRTY_VALUES,
            Column.EVENT_DIRTY,
        ], this.mapping.clone(), (this.mapping = mapping));
    };
    TestColumn.prototype.getColor = function (row) {
        return NumberColumn.prototype.getColor.call(this, row);
    };
    TestColumn.prototype.isFiltered = function () {
        return NumberColumn.prototype.isFiltered.call(this);
    };
    TestColumn.prototype.getFilter = function () {
        return NumberColumn.prototype.getFilter.call(this);
    };
    TestColumn.prototype.setFilter = function (value) {
        NumberColumn.prototype.setFilter.call(this, value);
    };
    // filter(row: IDataRow) {
    //   return NumberColumn.prototype.filter.call(this, row);
    // }
    /** @internal */
    TestColumn.prototype.isNumberIncluded = function (filter, value) {
        if (!filter) {
            return true;
        }
        if (Number.isNaN(value)) {
            return !filter.filterMissing;
        }
        return !((isFinite(filter.min) && value < filter.min) ||
            (isFinite(filter.max) && value > filter.max));
    };
    /**
     * filter the current row if any filter is set
     * @param row
     * @returns {boolean}
     */
    // TODO: customize filter: max, min, median...
    TestColumn.prototype.filter = function (row) {
        // currently it checks, if the median is within the range
        // const value = this.getRawNumber(row);
        var value = this.get_advanced_value(EAdvancedSortMethod.median, this.getRawValue(row)[0]["value"]);
        return this.isNumberIncluded(this.getFilter(), value);
    };
    TestColumn.prototype.clearFilter = function () {
        return NumberColumn.prototype.clearFilter.call(this);
    };
    var TestColumn_1;
    TestColumn.EVENT_MAPPING_CHANGED = NumberColumn.EVENT_MAPPING_CHANGED;
    TestColumn.EVENT_COLOR_MAPPING_CHANGED = NumberColumn.EVENT_COLOR_MAPPING_CHANGED;
    TestColumn.EVENT_SORTMETHOD_CHANGED = NumberColumn.EVENT_SORTMETHOD_CHANGED;
    TestColumn.EVENT_FILTER_CHANGED = NumberColumn.EVENT_FILTER_CHANGED;
    TestColumn = TestColumn_1 = __decorate([
        toolbar("rename", "filterNumber", "sort", "sortBy"),
        dialogAddons("sort", "sortNumbers"),
        SortByDefault("descending")
        //@ts-ignore
    ], TestColumn);
    return TestColumn;
}(MapColumn));
export { TestColumn };
//# sourceMappingURL=TestColumn.js.map