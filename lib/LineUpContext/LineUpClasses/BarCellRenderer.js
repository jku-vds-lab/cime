import { colorOf, DEFAULT_COLOR, ERenderMode, isNumberColumn, isNumbersColumn, renderMissingCanvas, renderMissingDOM, } from "lineupjs";
import { noRenderer, setText, CANVAS_HEIGHT, } from "./utils";
// https://github.com/lineupjs/lineupjs/blob/master/src/renderer/BarCellRenderer.ts
var BarCellRenderer = /** @class */ (function () {
    /**
     * flag to always render the value
     * @type {boolean}
     */
    function BarCellRenderer(renderValue) {
        if (renderValue === void 0) { renderValue = false; }
        this.renderValue = renderValue;
        this.title = "Bar";
    }
    BarCellRenderer.prototype.canRender = function (col, mode) {
        return (mode === ERenderMode.CELL && isNumberColumn(col) && !isNumbersColumn(col));
    };
    BarCellRenderer.prototype.create = function (col, context, imposer) {
        var width = context.colWidth(col);
        return {
            template: "<div title=\"\">\n          <div class=\"lu-bar-label\" style='background-color: " + DEFAULT_COLOR + "'>\n            <span " + (this.renderValue ? "" : "class=\"lu-hover-only\"") + "></span>\n          </div>\n        </div>",
            update: function (n, d) {
                var value = col.getNumber(d);
                var missing = renderMissingDOM(n, col, d);
                var w = isNaN(value) ? 0 : Math.round(value * 10000) / 100;
                var title = col.getLabel(d);
                n.title = title;
                var bar = n.firstElementChild;
                bar.style.width = missing ? "100%" : w + "%";
                var color = colorOf(col, d, imposer, value);
                //@ts-ignore
                bar.style.backgroundColor = missing ? null : color;
                setText(bar.firstElementChild, title);
                var item = bar.firstElementChild;
                setText(item, title);
                item.style.color = "black";
                // adaptDynamicColorToBgColor(item, color || DEFAULT_COLOR, title, w / 100);
            },
            render: function (ctx, d) {
                if (renderMissingCanvas(ctx, col, d, width)) {
                    return;
                }
                var value = col.getNumber(d);
                ctx.fillStyle = colorOf(col, d, imposer, value) || DEFAULT_COLOR;
                var w = width * value;
                ctx.fillRect(0, 0, isNaN(w) ? 0 : w, CANVAS_HEIGHT);
            },
        };
    };
    BarCellRenderer.prototype.createGroup = function () {
        return noRenderer;
    };
    BarCellRenderer.prototype.createSummary = function () {
        return noRenderer;
    };
    return BarCellRenderer;
}());
export default BarCellRenderer;
//# sourceMappingURL=BarCellRenderer.js.map